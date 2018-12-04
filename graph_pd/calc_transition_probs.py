# from igraph import Graph
import numpy as np
import multiprocessing as mp
from functools import partial
import os
import sys

from pd_graph_lib import run_death_birth_simulation_til_fixation, prisoners_dilemma_averaged

def get_fitness(cell_types,i):
    return 1+DELTA*game(cell_types[i],cell_types[neighbours[i]],b,c)

def read_neighbour_list(fname):
    neighbours = []
    with open(fname,'r') as f:
        for line in f:
            neighbours.append(np.fromstring(line,dtype=int,sep=' '))
    return neighbours

def calc_transition_probs(cell_types):
    fitnesses = [get_fitness(cell_types,i) for i in range(N)]
    T_m = 1./N*sum((1-cell_types[i])*cell_types[j]*fitnesses[i]/sum(fitnesses[k] for k in neighbours[j]) 
                for i in range(N)
                for j in neighbours[i])
    T_p = 1./N*sum((1-cell_types[j])*cell_types[i]*fitnesses[i]/sum(fitnesses[k] for k in neighbours[j]) 
                for i in range(N)
                for j in neighbours[i])
    return sum(cell_types),T_p,T_m

def run_sim(b,i):
    rand = np.random.RandomState()
    fix,history = run_death_birth_simulation_til_fixation(neighbours,DELTA,game,(b,1.0),timend,timestep,rand,return_fix = True,migration_strength=None)
    if fix == 1:
        data = [calc_transition_probs(cell_types)
                for cell_types in history]
        np.savetxt('%s/data_%d'%(outdir,i),data,fmt=('%3d','%.6e','%.6e'))

DELTA,game = 0.025,prisoners_dilemma_averaged        
c=1.0

nlist_file = 'hex_graph'
neighbours = read_neighbour_list(nlist_file)
neighbours = [np.array(n,dtype=int) for n in neighbours]

N = len(neighbours)
timend,timestep = 100000.,10.
sim_runs = int(sys.argv[1])
outdir = 'transition_probs/raw_data'
if not os.path.exists(outdir): # if the outdir doesn't exist create it
     os.makedirs(outdir)
b=8.
cpunum=mp.cpu_count()
pool = mp.Pool(processes=cpunum-1,maxtasksperchild=1000)
pool.map(partial(run_sim,b),range(sim_runs))
pool.close()
pool.join()

