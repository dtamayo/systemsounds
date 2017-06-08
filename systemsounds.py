from midiutil import MIDIFile
import rebound
from rebound import hash as rebhash
import numpy as np
import matplotlib.pyplot as plt
from subprocess import call
from itertools import repeat
import PIL # reminder that this is a requirement
from scipy.misc import imread
import warnings
warnings.filterwarnings("ignore")

def calc_midi_notes(particles, ref_note, ref_ID): # creates midi notes by scaling orbital frequencies using ref_note for particles[ref_ID] as the reference
    # 12 notes between octaves, freq = f0*2**(n/12). 
    # n = 12 log_2(freq/f0)
    # star, then planets from inside out
    midinotes = [0] # placeholder for star
    for p in particles[1:]:
        midinote = ref_note+12*np.log(particles[ref_ID].P/p.P)/np.log(2)
        midinotes.append(int(np.round(midinote)))
    return midinotes  

def rescale_time(sim, timescale): # Rescale time by timescale, i.e. timescale time units will now correspond to 1 time units
    ps = sim.particles
    sim.G *= timescale**2
    sim.dt /= timescale
    for p in ps:
        p.vx *= timescale
        p.vy *= timescale
        p.vz *= timescale

def prepend_to_heartbeat(sim, func):
    try:    # can't store sim.heartbeat b/c ctypes returns new object each time so we'd always get latest hb
        oldhb = sim._pyheartbeat                      
    except AttributeError:
        oldhb = lambda s: None
    def heartbeat(reb_sim):
        oldhb(reb_sim)
        func(reb_sim)
    sim._pyheartbeat = heartbeat                            # store python function in sim
    sim.heartbeat = heartbeat                               # update ctypes function wrapper in sim

def copysim(sim):                                     # should eventually add better version to REBOUND
    sim2 = rebound.Simulation()
    sim2.G = sim.G
    sim2.t = sim.t
    for p in sim.particles:
        sim2.add(p)
    return sim2

class EventRecorder(object):
    def __init__(self, sim, rootfunc, targets=[None]):
        self.events = []
        self._oldvals = {}
        self.targets = targets
        self.add_event_recorder_to_heartbeat(sim, rootfunc)
    def add_event_recorder_to_heartbeat(self, sim, rootfunc):
        def check_for_root_crossings(reb_sim):
            sim = reb_sim.contents
            ps = sim.particles
            for target in self.targets:
                val = rootfunc(sim, target)
                if self._oldvals[target] is not None:           # not first call
                    if self._oldvals[target] < 0 and val > 0:   # crossed from negative to positive
                        sim_root_crossing = self.bisection(sim, rootfunc, target)
                        self.process_event(sim_root_crossing, target)
                self._oldvals[target] = val
        prepend_to_heartbeat(sim, check_for_root_crossings)

    def process_event(self, sim, target, params={}):
        params['t'] = sim.t
        params['target'] = target
        self.events.append(params)

    def bisection(self, sim, rootfunc, target, epsilon=1.e-6): # bisection to find crossing time
        sim2 = copysim(sim)
        oldt = sim.t # need to go back from overshot t to previous value
        newt = sim.t - sim.dt_last_done/2. 
        sim2.dt *= -1
        oldval = rootfunc(sim2, target) 
        while (abs(newt - oldt)/oldt > epsilon):
            midt = (newt + oldt)/2.
            sim2.integrate(midt)
            val = rootfunc(sim2, target)
            if oldval*val < 0.: # switched sign
                newt = oldt # go back to prev value
                sim2.dt *= -0.3
            else: # keep integrating toward newt
                sim2.dt *= 0.3
            oldt = midt # next iteration starts at midt
            oldval = val
        return sim2

    @property
    def targets(self):
        return self._targets
    @targets.setter
    def targets(self, iterator):
        self._targets = iterator
        for target in iterator:
            if target not in self._oldvals.keys():
                self._oldvals[target] = None

class FrameRecorder(EventRecorder):
    def __init__(self, sim, time_per_sec, fps=30):
        try:
            call("rm -f ./tmp/*", shell=True)
        except:
            pass
        self.fps = fps
        self.time_per_sec = time_per_sec
        self._frame_ctr = 0
        self._last_frame_time = sim.t
        
        def root_func(sim, target=None):
            return sim.t - self._last_frame_time - 1./self.fps
        super(FrameRecorder, self).__init__(sim, root_func)

    def process_event(self, frame_sim, target=None):
        filename = "tmp/binaries/frame"+str(self._frame_ctr)+".bin"
        frame_sim.save(filename)
        self._last_frame_time = frame_sim.t
        
        params={}
        for key in self.__dict__.keys():
            if key is not "events":
                params[key] = self.__dict__[key]
        self._frame_ctr += 1
                            
        super(FrameRecorder, self).process_event(frame_sim, target, params)
            
def write_png(params):
    fig_ctr, time, filename, time_per_beat, color, showparticles, showtransits, showconjunctions, conjunctions, background, transparent = params
    coloriterator = [color[i] for i in showparticles]
    sim = rebound.Simulation.from_file(filename)
    sim.t=0
    rescale_time(sim, time_per_beat)
    sim.integrate(time)
    ps = sim.particles
    
    lw=3
    fadetimescale = sim.particles[-1].P/3. # for conjunctions
    refsize=25*lw # this is what REBOUND uses for size of circles in call to plt.scatter

    fig = rebound.OrbitPlot(sim, figsize=(8,8), color=coloriterator, lw=lw, plotparticles=showparticles)
    ax = fig.axes[0]
    ax.axis('off')
        
    for i in showparticles:
        p = ps[i]
        ax.scatter(p.x, p.y, s=refsize, color=color[i], marker='o', zorder=4)
    for i in showtransits:
        p = ps[i]
        scale=p.a/3 # length scale for making dots bigger
        if p.x > 0 and np.abs(p.y)/scale < 1:
            ax.scatter(p.x, p.y, s=refsize*(1+6*np.exp(-np.abs(p.y)/scale)),color=color[i], marker='o', zorder=5)
    
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    cscale = 10*xlim[1]
    if showconjunctions:
        nearby_conjunctions = [conjunction for conjunction in conjunctions if time - conjunction[0] < fadetimescale]
        for conjunction in nearby_conjunctions:
            j, x, y = conjunction[1], conjunction[2], conjunction[3]
            if j in showconjunctions and j+1 in showconjunctions:
                ax.plot([0, cscale*x], [0,cscale*y], lw=5, color=color[j], alpha=max(1.-(time-conjunction[0])/fadetimescale,0.), zorder=1)
       
    if background:
        bkg = imread('images/US_background_image.png')
        ax.imshow(bkg, zorder=0, extent=xlim+ylim)
    fig.savefig('tmp/pngs/{0:0=5d}.png'.format(fig_ctr), transparent=transparent,dpi=300)
    plt.close(fig)  

class Event():
    def __init__(self, t, target):
        self.t = t
        self.target = target

class System(rebound.Simulation):
    def __init__(self, fps=30):
        super(System, self).__init__()
        self.initialize(fps=30)

    @classmethod
    def from_file(cls, filename):
        sim = rebound.Simulation.from_file(filename)
        sim.__class__ = cls
        sim.initialize(fps=30)
        return sim

    def initialize(self, fps=30):
        try:
            call("rm -f ./tmp/*", shell=True)
        except:
            pass
        self.t = 0
        self.fps = fps
        self._frame_ctr = 0
        self._fig_timer = 0.
        self.time_elapsed = 0
        
        self.fig_params = []
        self.conjunctions = []
        self.transits = []
        self.frames = []
        self.tempo = []

        self.recordtransits = []
        self.recordconjunctions = []
        
        self.time_per_sec = None
        self.set_heartbeat()
    def set_heartbeat(self):
        def heartbeat_wrapper(self):
            def heartbeat(reb_sim):
                time_elapsed = self.dt_last_done/self.time_per_sec 
                self.time_elapsed += time_elapsed 
                self._fig_timer += time_elapsed
                if self._fig_timer > 1./self.fps:
                    filename = "tmp/frame"+str(self._frame_ctr)+".bin"
                    self.save(filename)
                    self._frame_ctr += 1
                    self.frames.append(Frame(self.t, filename))
                    #print(self._frame_ctr, self._frame_ctr/self.fps, self.t/self.time_per_sec)
                    self._fig_timer -= 1./self.fps
            return heartbeat
        hb_wrapper = heartbeat_wrapper(self)
        self.heartbeat = hb_wrapper
        self._hbeat = hb_wrapper 

        #add_rootfinder_to_heartbeat(self, lambda sim, i: sim.particles[i].y, lambda s: s.transits, lambda s: s.recordtransits)
        #add_rootfinder_to_heartbeat(self, lambda sim, i: np.sin(sim.particles[i].theta - sim.particles[i+1].theta), lambda s: s.conjunctions, lambda s: s.recordconjunctions)
   
    @property
    def time_per_sec(self): # time_per_sec is beats per second, so just bpm/60
        return self._time_per_sec
    @time_per_sec.setter
    def time_per_sec(self, value):
        self._time_per_sec = value
    def write_images(self):
        call("rm -f tmp/pngs/*", shell=true)
        pool = rebound.InterruptiblePool()
        for a in self.fig_params:
            a.append(None)
            a.append(False)
        res = pool.map(write_png, self.fig_params)
    def write_movie(self, moviename, midiname=None):
        fps = 30
        try:
            call("rm -f {0}.mp4".format(moviename), shell=True)
        except:
            pass
        
        if midiname:
            try:
                call("rm -f ./tmp/{0}.wav".format(midiname), shell=True)
                call("rm -f ./tmp/{0}cut.wav".format(midiname), shell=True)
            except:
                pass
            call("timidity -Ow ./{0}.mid -o ./tmp/{0}.wav --preserve-silence".format(midiname), shell=True)
            call("ffmpeg -t {0} -i ./tmp/{1}.wav ./tmp/{1}cut.wav".format(self.time_elapsed, midiname), shell=True)
            call("ffmpeg -r {0} -i ./tmp/pngs/%05d.png -i tmp/{1}cut.wav -c:v libx264 -pix_fmt yuv420p -c:a libfdk_aac -b:a 192k -shortest {2}.mp4".format(fps, midiname, moviename), shell=True)   
        else:
            try:
                call("rm -f {0}".format(moviename), shell=True)
            except:
                pass
            call("ffmpeg -r {0} -i tmp/pngs/%05d.png -c:v libx264 -pix_fmt yuv420p {1}.mp4".format(fps, moviename), shell=True)        
