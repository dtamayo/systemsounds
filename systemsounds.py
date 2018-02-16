from midiutil import MIDIFile
import rebound
import numpy as np
import matplotlib.pyplot as plt
from subprocess import call
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
        oldhb = sim._hb
    except AttributeError:
        oldhb = lambda s: None
    def heartbeat(reb_sim):
        oldhb(reb_sim)
        func(reb_sim)
    sim.heartbeat = heartbeat                               # update ctypes function wrapper in sim

def copysim(sim):                                           # should eventually add better version to REBOUND
    sim2 = rebound.Simulation()
    sim2.G = sim.G
    sim2.t = sim.t
    for p in sim.particles:
        sim2.add(p)
    return sim2

class EventRecorder(object):
    def __init__(self, sim, rootfunc, targets=None, verbose=False):
        self.events = []
        self._oldvals = {}
        self.targets = range(1, sim.N) if targets is None else targets
        self.add_event_recorder_to_heartbeat(sim, rootfunc)
        self.verbose = verbose
    def add_event_recorder_to_heartbeat(self, sim, rootfunc):
        def check_for_root_crossings(reb_sim):
            sim = reb_sim.contents
            ps = sim.particles
            for target in self.targets:
                val = rootfunc(sim, target)
                #print(sim.t, self._oldvals[None], val, type(target))
                if self._oldvals[target] is not None and self._oldvals[target] < 0 and val >= 0:   # not first call, and crossed from negative to positive
                    sim_root_crossing = self.bisection(sim, self._oldvals[target], rootfunc, target)
                    self.process_event(sim_root_crossing, target)
                    self._oldvals[target] = -1. # set oldval to negative number so we get a root crossing if another event happens within next timestep
                    if self.verbose:
                        print("{0} event at t = {1}".format(self.__class__.__name__, sim_root_crossing.t))
                else:
                    self._oldvals[target] = rootfunc(sim, target)
        prepend_to_heartbeat(sim, check_for_root_crossings)

    def process_event(self, sim, target):
        params={'time':sim.t, 'target':target}
        for key in self.__dict__.keys():
            if key.startswith('_') is False and key is not "events":
                params[key] = self.__dict__[key]
        self.events.append(params)

    def bisection(self, sim, val, rootfunc, target, epsilon=1.e-6): # bisection to find crossing time
        sim2 = copysim(sim)
        oldt = sim.t # need to go back from overshot t to previous value
        newt = sim.t - sim.dt_last_done
        sim2.dt *= -1
        oldval = rootfunc(sim2, target) 
        if abs(oldval/val) > epsilon: # check edge case where oldval is at 0 initially (use val for scale)
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
    def __init__(self, sim, time_per_sec, fps=30, plotparticles=None, verbose=False):
        try:
            call("rm -f ./tmp/*", shell=True)
        except:
            pass
        self.fps = fps
        self.time_per_sec = time_per_sec
        self.frame_ctr = 0
        self.elapsed_time = 0. # movie time in seconds
        self._last_frame_time = sim.t
        self.plotparticles = range(1, sim.N) if plotparticles is None else plotparticles

        def update_elapsed_time(reb_sim):
            sim = reb_sim.contents
            self.elapsed_time += sim.dt_last_done/self.time_per_sec
            if sim.dt > self.time_per_sec/self.fps: # check timestep is not longer than time between frames
                sim.dt = self.time_per_sec/self.fps/np.sqrt(2.) # make timestep shorter than time between frames (random irrational close to 1)
                sim.ri_ias15.epsilon = 0       # set constant timestep in ias15 so doesn't increase
        prepend_to_heartbeat(sim, update_elapsed_time)

        def root_func(sim, target=None):
            return sim.t - self._last_frame_time - self.time_per_sec/self.fps
        super(FrameRecorder, self).__init__(sim, root_func, targets=[None], verbose=verbose) # no individual targets for timer, so pass iterator with single entry

    def process_event(self, frame_sim, target=None):
        self.filename = "tmp/binaries/frame"+str(self.frame_ctr)+".bin"
        frame_sim.save(self.filename)
        self._last_frame_time = frame_sim.t
        super(FrameRecorder, self).process_event(frame_sim, target)
        self.frame_ctr += 1
