import pandas as pd
import rebound 
import numpy as np

def add_planet(sim, planet, T1):
    P = planet['pl_orbper']/365.25
    e = planet['pl_orbeccen']
    pomega = planet['pl_orbeccen']/(2.*np.pi)
    t_transit = (planet['pl_tranmid']-T1)/365.25 # time of transit relative to first planet, defined as t=0
    for el in [P, e, pomega, t_transit]:
        if np.isnan(el) is True:
            el = 0.
    sim.integrate(t_transit)
    sim.add(P=P, e=e, pomega=pomega, theta=0) # add new planet at true longitude=0 (transiting)
    
def make_sim(df, name):
    system = df[df['pl_hostname'] == name]
    if system.shape[0] == 0:
        raise AttributeError("System {0} not found in NASA exoplanet archive".format(name))
    
    p1 = system.iloc[0]
    starmass = p1['st_mass']
    sim = rebound.Simulation()
    sim.G = 4*np.pi**2
    sim.add(m=starmass)
    
    T1 = p1['pl_tranmid'] # time of mid-transit for first planet
    for i, planet in system.iterrows():
        add_planet(sim, planet, T1)

    Ps = [sim.particles[i].P for i in range(1, sim.N)]
    Pmin = np.array(Ps).min()
    sim.integrator="whfast"
    sim.dt = Pmin*0.07
    return sim
