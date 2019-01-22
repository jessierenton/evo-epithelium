from multiprocessing import Pool  #parallel processing
import multiprocessing as mp
import structure
from structure.global_constants import *
from structure.cell import Tissue, BasicSpringForceNoGrowth
import structure.initialisation as init
import sys
import os
import numpy as np
import libs.pd_lib_neutral as lib
import libs.data as data


def calc_interactions(tissue,mutant_index,n):
    """treats all cells with ancestor 'mutant_index' as cooperators
    returns:
        n (int): size of clone
        I_CC/I_CD (ints): number of cooperator-cooperator/defector interactions in population
        W_CC/W_CD (floats): number of cooperator-cooperator/defector interactions in pop. weighted by neighbour number    
    """
    neighbours = tissue.mesh.neighbours
    types = tissue.properties['ancestor']==mutant_index
    I_CC,I_CD,W_CC,W_CD,N_D = 0,0,0.,0.,0
    for ctype,cell_neighbours in zip(types,neighbours):
        if ctype:
            Cneigh,neigh = float(sum(types[cell_neighbours])),float(len(cell_neighbours))
            I_CC += Cneigh
            I_CD += neigh - Cneigh
            W_CC += Cneigh/neigh
            W_CD += (neigh-Cneigh)/neigh
    return [n,I_CC,I_CD,W_CC,W_CD]

def run_sim(i):
    """run a single simulation and save interaction data for each clone"""
    rand = np.random.RandomState()
    dt=0.005*-50./MU
    tissue = lib.initialise_tissue_ancestors(l,dt,10.,10.,rand,MU)
    tissue.properties['ancestor']=np.arange(l*l)
    if init_timend is not None: tissue = lib.run(tissue,lib.simulation_ancestor_tracking(tissue,dt,init_timend/dt,init_timend/dt,rand),init_timend/dt,init_timend/dt)[-1]
    data = [calc_interactions(tissue,mutant_index,n)
                for tissue in lib.run_generator(lib.simulation_ancestor_tracking(tissue,dt,timend/dt,timestep/dt,rand,til_fix=True),timend/dt,timestep/dt)
                for mutant_index,n in enumerate(np.bincount(tissue.properties['ancestor'])) if n>=n_min]                
    np.savetxt('%s/data_%d'%(outdir,i),data,fmt=('%4d','%4d','%4d','%4.6f','%4.6f'))
    return None

l = 10 # population size N = l*l
init_timend = 10. # initial simulation time to equilibrate 
timestep = 12. # timesteps at which to calc interaction data (hours)
timend = 10000. # length of simulation (hours)
sim_runs = int(sys.argv[1]) # number of sims to run taken as command line arg
MU = float(sys.argv[2]) #spring constant
n_min = 1 
outdir = 'interaction_data/supp_vary_MU/MU%d/raw_data'%MU
if not os.path.exists(outdir): # if the outdir doesn't exist create it 
     os.makedirs(outdir)


# run simulations in parallel 
cpunum=mp.cpu_count()
pool = Pool(processes=cpunum-1,maxtasksperchild=1000)
pool.map(run_sim,range(sim_runs))
pool.close()
pool.join()
