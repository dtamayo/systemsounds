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

"""
integrate:
    store events
    recordtransits
    transits = []
    frames = []
    Event:
        init(t)
    Transit
        init(particle)
    Conjunction
        init(pair)
    self.midi.addNote(track, ps[j].index, self.notes[j], t, duration, self.velocities[j])
    self.conjunctions.append((self.sim.t, j, ps[j].x, ps[j].y))
    self.midi.addNote(track, N, self.conjunction_notes[j], t, duration, self.conjunction_velocities[j]) # add to track above all planets
    hash?
"""


def rescale_time(sim, timescale): # Rescale time by timescale, i.e. timescale time units will now correspond to 1 time units
    ps = sim.particles
    sim.G *= timescale**2
    sim.dt /= timescale
    for p in ps:
        p.vx *= timescale
        p.vy *= timescale
        p.vz *= timescale

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
    def __init__(self, filename):
        self.filename = filename

class Frame(Event):
    def __init__(self, filename, fig_ctr):
        Event.__init__(self, filename)
        self.fig_ctr = fig_ctr

class Transit(Event):
    def __init__(self, filename, particleID):
        Event.__init__(self, filename)
        self.particleID = particleID

class Conjunction(Event):
    def __init__(self, filename, innerID, outerID):
        Event.__init__(self, filename)
        self.innerID = innerID
        self.outerID = outerID

class System(rebound.Simulation):
    def __init__(self, time_per_sec=1, dt=None, dt_epsilon=1.e-5, outer_midi_note=48, fps=30, exact_midi_times=True):
        super(System, self).__init__()
        self.initialize(time_per_sec=1, dt=None, dt_epsilon=1.e-5, outer_midi_note=48, fps=30, exact_midi_times=True)

    @classmethod
    def from_file(cls, filename):
        sim = rebound.Simulation.from_file(filename)
        sim.__class__ = cls
        sim.initialize(time_per_sec=1, dt=None, dt_epsilon=1.e-5, outer_midi_note=48, fps=30, exact_midi_times=True)
        return sim

    def initialize(self, time_per_sec=1, dt=None, dt_epsilon=1.e-5, outer_midi_note=48, fps=30, exact_midi_times=True):
        try:
            call("rm -f ./tmp/*", shell=True)
        except:
            pass
        self.midi = MIDIFile(adjust_origin=True) # One track, defaults to format 1 (tempo track automatically created)
        #self.sim = sim
        self.exact_midi_times = exact_midi_times
        self.t = 0
        self.dt_epsilon = dt_epsilon
        
        self.time_per_sec = time_per_sec
        self.fps = fps
        self.fig_ctr = 0
        self.event_ctr = 0
        self.time_elapsed = 0
        self.time_per_fig = 1./self.fps
        
        #self.notes = self.calc_midi_notes(outer_midi_note)
        #self.velocities = [100 for i in range(self.sim.N)]
        #self.conjunction_notes = [12 for i in range(self.sim.N)] # C1
        #self.conjunction_velocities = [100 for i in range(self.sim.N)]
        
        #if not dt:
        #    self.dt = self.sim.particles[1].P/5.
        #else:
        #    self.dt = dt
        self.fig_params = []
        self.conjunctions = []
        self.transits = []
        self.frames = []

        self._yprev = [None]*self.N
    
    def calc_midi_notes(self, outer_midi_note):
        # 12 notes between octaves, freq = f0*2**(n/12). 
        # n = 12 log_2(freq/f0)
        # star, then planets from inside out
        ps = self.sim.particles
        midinotes = [0] # placeholder for star
        for p in ps[1:]:
            midinote = outer_midi_note+12*np.log(ps[-1].P/p.P)/np.log(2)
            midinotes.append(int(np.round(midinote)))
        return midinotes  
    def make_tuple(self, arg):
        N = self.sim.N
        if arg == True:
            return tuple(range(1,N))
        elif arg == False:
            return ()
        else:
            return tuple(arg)
    def copysim(self):
        sim = rebound.Simulation()
        sim.G = self.sim.G
        sim.t = self.sim.t
        for p in self.sim.particles:
            sim.add(p)
        return sim

    def find_exact_crossing_time(self, get_val, prevt): # do bisection to find exact crossing time
        sim2 = self.copysim()
        oldt = self.sim.t # need to go back from overshot t to previous value
        newt = prevt 
        sim2.dt *= -1
        oldval = get_val(sim2) #ps[j].y
        while (abs(newt - oldt)> self.dt_epsilon):
            midt = (newt + oldt)/2.
            #print(oldt, newt, midt, get_val(sim2))
            sim2.integrate(midt)
            if oldval*get_val(sim2) < 0.: # switched sign
                newt = oldt # go back to prev value
                sim2.dt *= -0.3
            else: # keep integrating toward newt
                sim2.dt *= 0.3
            oldt = midt # next iteration starts at midt
            oldval = get_val(sim2)
        return sim2

    def findtransits(self, yprev, prevt):
        ps = self.sim.particles
        for i in range(1, self.sim.N):
            if yprev[i] < 0 and ps[i].y > 0: # Crossed x axis
                transitsim = self.find_exact_crossing_time(lambda sim: sim.particles[i].y, prevt)
                print("Transit: {0}\t{1}".format(transitsim.t, i))
                #self.midi.addNote(track, ps[pid].index, self.notes[i], t, duration, self.velocities[i])
                filename = ''
                #filename = "tmp/event"+str(self.event_ctr)+".bin"
                #self.sim.save(filename)
                #self.event_ctr += 1
                self.transits.append(Transit(filename, i))
            yprev[i] = ps[i].y
        return yprev

    def findconjunctions(self, sinthetaprev, prevt):
        ps = self.sim.particles
        for i, pair in enumerate(self.recordconjunctions):
            innerID, outerID = pair
            sintheta = np.sin(ps[outerID].theta-ps[innerID].theta)
            if sinthetaprev[i] > 0 and sintheta < 0:
                if self.exact_midi_times:
                    t = self.find_exact_crossing_time(lambda sim: np.sin(sim.particles[outerID].theta-sim.particles[innerID].theta), prevt)
                else:
                    t = self.sim.t
                print('Conjunction: {0}\t{1}'.format(t, innerID))
                #self.conjunctions.append((self.sim.t, innerID, ps[j].x, ps[j].y))
                filename = "tmp/event"+str(self.event_ctr)+".bin"
                self.sim.save(filename)
                self.event_ctr += 1
                self.conjunctions.append(Conjunction(filename, innerID, outerID))
                #self.midi.addNote(track, N, self.conjunction_notes[innerID], t, duration, self.conjunction_velocities[innerID]) # add to track above all planets
            sinthetaprev[i] = sintheta
        return sinthetaprev
    
    def copysim2(sim):
        sim2 = rebound.Simulation()
        sim2.G = sim.G
        sim2.t = sim.t
        for p in sim.particles:
            sim2.add(p)
        return sim2
    
    def find_exact_crossing_time2(sim, f): # do bisection to find exact crossing time
        sim2 = copysim2(sim)
        oldt = self.sim.t # need to go back from overshot t to previous value
        newt = prevt 
        sim2.dt *= -1
        oldval = get_val(sim2) #ps[j].y
        while (abs(newt - oldt)> self.dt_epsilon):
            midt = (newt + oldt)/2.
            #print(oldt, newt, midt, get_val(sim2))
            sim2.integrate(midt)
            if oldval*get_val(sim2) < 0.: # switched sign
                newt = oldt # go back to prev value
                sim2.dt *= -0.3
            else: # keep integrating toward newt
                sim2.dt *= 0.3
            oldt = midt # next iteration starts at midt
            oldval = get_val(sim2)
        return sim2

    def findtransits2(sim):
        ps = sim.particles
        for i in range(1, sim.N):
            if sim._yprev[i] < 0 and ps[i].y > 0: # Crossed x axis
                transitsim = find_exact_crossing_time2(sim, lambda s: s.particles[i].y)
                print("Transit: {0}\t{1}".format(transitsim.t, i))
                #self.midi.addNote(track, ps[pid].index, self.notes[i], t, duration, self.velocities[i])
                filename = ''
                #filename = "tmp/event"+str(self.event_ctr)+".bin"
                #self.sim.save(filename)
                #self.event_ctr += 1
                self.transits.append(Transit(filename, i))
            sim._yprev[i] = ps[i].y

    def hb(sim):
        findtransits2(sim)

    '''
    def integrate(self, tmax, color=True, duration=1, track=0, planetentrance=False):
        new_entrance=False  #if planetentrance==True this will only show innermost planet once it transits for the first time
        
        N=self.sim.N
        ps = self.sim.particles
        yprev = np.zeros(self.sim.N)
        #sinthetaprev = np.zeros(len(self.recordconjunctions))
        while self.sim.t < tmax:
            prevt = self.sim.t
            self.sim.integrate(self.sim.t+self.dt)
            #print(self.sim.t, self.sim.particles[1].x, self.sim.particles[1].y)
            self.time_elapsed += self.dt/self.time_per_sec #self.dt/self.bpm*60.
            
            yprev = self.findtransits(yprev, prevt) 
            #sinthetaprev = self.findconjunctions(sinthetaprev, prevt)

            if self.time_elapsed/self.time_per_fig > self.fig_ctr + 1:
                if len(self.recordtransits)>1 and planetentrance==True and new_entrance==False:
                    #print('Waiting for new entrance:', self.sim.t,tuple(self.showparticles))
                    print('Figure waiting for new entrance: {0}\t{1}'.format(self.fig_ctr, self.sim.t))
                    filename = "tmp/event"+str(self.event_ctr)+".bin"
                    self.sim.save(filename)
                    self.event_ctr += 1
                    self.frames.append(Frame(filename, self.fig_ctr))
                    #self.fig_params.append([self.fig_ctr, self.sim.t, color, tuple(self.showparticles)[1:], tuple(self.showtransits)[1:], tuple(self.showconjunctions), tuple(self.conjunctions)])
                else:
                    print('Figure: {0}\t{1}'.format(self.fig_ctr, self.sim.t))
                    #self.fig_params.append([self.fig_ctr, self.sim.t, color, tuple(self.showparticles), tuple(self.showtransits), tuple(self.showconjunctions), tuple(self.conjunctions)])
                    filename = "tmp/event"+str(self.event_ctr)+".bin"
                    self.sim.save(filename)
                    self.frames.append(Frame(filename, self.fig_ctr))
                self.fig_ctr += 1
    '''
    def change_tempo(self, bpm):
        self.bpm = bpm
        self.time_per_sec = self.bpm/60.
        self.midi.addTempo(0, self.sim.t, self.bpm) 
    def write_midi(self, midiname):
        with open("./"+midiname+".mid", "wb") as f:
            self.midi.writeFile(f)
    def write_images(self):
        call("rm -f tmp/pngs/*", shell=True)
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
        
