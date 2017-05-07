from midiutil import MIDIFile
import rebound
import numpy as np
import matplotlib.pyplot as plt
from subprocess import call
from itertools import repeat
import PIL # reminder that this is a requirement
from scipy.misc import imread
import warnings
warnings.filterwarnings("ignore")

def set_time_per_beat(sim, time_per_beat): # makes sim.t run in units of the outer planet orbit = one beat
    ps = sim.particles
    sim.G = time_per_beat**2
    sim.dt /= time_per_beat
    for p in ps:
        p.vx *= time_per_beat
        p.vy *= time_per_beat
        p.vz *= time_per_beat

def write_png(params):
    fig_ctr, time, filename, time_per_beat, color, showplanets, showtransits, showconjunctions, conjunctions, background, transparent = params
    coloriterator = [color[i] for i in showplanets]
    sim = rebound.Simulation.from_file(filename)
    sim.t=0
    set_time_per_beat(sim, time_per_beat)
    sim.integrate(time)
    ps = sim.particles
    
    lw=3
    fadetimescale = sim.particles[-1].P/3. # for conjunctions
    refsize=25*lw # this is what REBOUND uses for size of circles in call to plt.scatter

    fig = rebound.OrbitPlot(sim, figsize=(8,8), color=coloriterator, lw=lw, plotparticles=showplanets)
    ax = fig.axes[0]
    ax.axis('off')
        
    for i in showplanets:
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

class System():
    def __init__(self, filename, bpm, time_per_beat=None, dt=None, dt_epsilon=1.e-5, outer_midi_note=48, fps=30, exact_midi_times=True):
        try:
            call("rm -f ./tmp/*", shell=True)
        except:
            pass
        self.midi = MIDIFile(adjust_origin=True) # One track, defaults to format 1 (tempo track automatically created)
        self.filename = filename
        self.sim = rebound.Simulation.from_file(filename)
        self.exact_midi_times = exact_midi_times
        if self.exact_midi_times:
            self.sim.integrator = "ias15"
        self.sim.t = 0
        self.dt_epsilon = dt_epsilon
        
        self.bpm = bpm
        self.fps = fps
        self.fig_ctr = 0
        self.time_elapsed = 0
        self.time_per_fig = 1./self.fps
        
        self.notes = self.calc_midi_notes(outer_midi_note)
        self.velocities = [100 for i in range(self.sim.N)]
        self.conjunction_notes = [12 for i in range(self.sim.N)] # C1
        self.conjunction_velocities = [100 for i in range(self.sim.N)]
        
        if time_per_beat:
            self.time_per_beat = time_per_beat
        else:
            self.time_per_beat = self.sim.particles[-1].P
        set_time_per_beat(self.sim, self.time_per_beat)
        self.change_tempo(bpm)
        if not dt:
            self.dt = self.sim.particles[1].P/5.
        else:
            self.dt = dt
        self.fig_params = []
        self.conjunctions = []
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
        return sim2.t

    def integrate(self, tmax, color=True, duration=1, track=0, playtransits=True, playconjunctions=False, showplanets=True, showtransits=True, showconjunctions=True, planetentrance=False):
        playtransits = self.make_tuple(playtransits)
        playconjunctions = self.make_tuple(playconjunctions)
        showplanets = self.make_tuple(showplanets)
        showtransits = self.make_tuple(showtransits)
        showconjunctions = self.make_tuple(showconjunctions)
        new_entrance=False  #if planetentrance==True this will only show innermost planet once it transits for the first time
        
        N=self.sim.N
        ps = self.sim.particles
        yprev = np.zeros(N)
        sinthetaprev = np.zeros(N)
        while self.sim.t < tmax:
            prevt = self.sim.t
            self.sim.integrate(self.sim.t+self.dt)
            #print(self.sim.t, self.sim.particles[1].x, self.sim.particles[1].y)
            self.time_elapsed += self.dt/self.bpm*60.
            for j in playtransits:
                if yprev[j] < 0 and ps[j].y > 0: # Crossed x axis
                    #if j==1:
                    #    print(self.sim.t)
                    if self.exact_midi_times:
                        t = self.find_exact_crossing_time(lambda sim: sim.particles[j].y, prevt)
                    else:
                        t = self.sim.t
                    self.midi.addNote(track, ps[j].index, self.notes[j], t, duration, self.velocities[j])
                    if j==playtransits[0]:  
                        new_entrance=True
                yprev[j] = ps[j].y
       
            for j in playconjunctions:
                if j+1 in playconjunctions:
                    sintheta = np.sin(ps[j+1].theta-ps[j].theta)
                    if sinthetaprev[j] > 0 and sintheta < 0:
                        if self.exact_midi_times:
                            t = self.find_exact_crossing_time(lambda sim: np.sin(sim.particles[j+1].theta-sim.particles[j].theta), prevt)
                        else:
                            t = self.sim.t
                        #print('conjunction', t, j)
                        self.conjunctions.append((self.sim.t, j, ps[j].x, ps[j].y))
                        self.midi.addNote(track, N, self.conjunction_notes[j], t, duration, self.conjunction_velocities[j]) # add to track above all planets
                    sinthetaprev[j] = sintheta
            if self.time_elapsed/self.time_per_fig > self.fig_ctr + 1:
                if len(playtransits)>1 and planetentrance==True and new_entrance==False:
                    #print('Waiting for new entrance:', self.sim.t,tuple(showplanets))
                    self.fig_params.append([self.fig_ctr, self.sim.t, self.filename, self.time_per_beat, color, tuple(showplanets)[1:], tuple(showtransits)[1:], tuple(showconjunctions), tuple(self.conjunctions)])
                else:
                    self.fig_params.append([self.fig_ctr, self.sim.t, self.filename, self.time_per_beat, color, tuple(showplanets), tuple(showtransits), tuple(showconjunctions), tuple(self.conjunctions)])
                self.fig_ctr += 1
    def change_tempo(self, bpm):
        self.bpm = bpm
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
        
