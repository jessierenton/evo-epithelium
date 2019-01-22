from multiprocessing import Pool
import numpy as np
import libs.pd_lib_neutral as lib #library for simulation routines
import libs.data as data
# import libs.plot as vplt #plotting library
from structure.global_constants import *
import structure.initialisation as init
from structure.cell import Tissue, BasicSpringForceNoGrowth

"""run a single voronoi tessellation model simulation"""

def run_sim(MU):
    rand = np.random.RandomState()
    dt=min(0.005,0.005*-50./MU)
    return lib.run_simulation_vary_MU(simulation,l,timestep,timend,rand,MU,dt)
    

l = 10 # population size N=l*l
timend = 500 # simulation time (hours)
timestep = 1. # time intervals to save simulation history

rand = np.random.RandomState()
simulation = lib.simulation_ancestor_tracking  #simulation routine imported from lib

MU_vals = -np.array((0.1,1,10,25,50,100,250),dtype=float)
pool = Pool(maxtasksperchild=1000)
for i in range(3):
    for MU,history in zip(MU_vals,pool.imap(run_sim,MU_vals)):
        data.save_force(history,'vary_MU_data/MU%.1f'%MU,index=i)
        data.save_neighbour_distr(history,'vary_MU_data/MU%.1f'%MU,index=i)

