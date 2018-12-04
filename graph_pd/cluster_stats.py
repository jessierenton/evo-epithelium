# from igraph import Graph
import numpy as np
import multiprocessing as mp
from functools import partial
import os
import sys

from pd_graph_lib import run_death_birth_simulation_neutral_clones


def read_neighbour_list(fname):
    neighbours = []
    with open(fname,'r') as f:
        for line in f:
            neighbours.append(np.fromstring(line,dtype=int,sep=' '))
    return neighbours

def calc_interactions(neighbours,cell_types,mutant_index,n):
    """treats all cells with ancestor 'mutant_index' as cooperators
    returns:
        n (int): size of clone
        I_CC/I_CD (ints): number of cooperator-cooperator/defector interactions in population
        W_CC/W_CD (floats): number of cooperator-cooperator/defector interactions in pop. weighted by neighbour number    
    """
    types = cell_types==mutant_index
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
    history = run_death_birth_simulation_neutral_clones(neighbours,timend,timestep,rand)
    data = [calc_interactions(neighbours,cell_types,mutant_index,n)
                for cell_types in history
                for mutant_index,n in enumerate(np.bincount(cell_types)) if n>0]
    np.savetxt('%s/data_%d'%(outdir,i),data,fmt=('%4d','%4d','%4d','%4.6f','%4.6f'))
            

nlist_file = 'hex_graph' #file from which to read graph structure
neighbours = read_neighbour_list(nlist_file)
neighbours = [np.array(n,dtype=int) for n in neighbours]

N = len(neighbours) #population size
timend,timestep = 100000.,20. #length of simulation (hours), timesteps at which to calc interaction data
sim_runs = int(sys.argv[1]) # number of sims to run taken as command line arg
outdir = 'interaction_data'
if not os.path.exists(outdir): # if the outdir doesn't exist create it
     os.makedirs(outdir)

cpunum=mp.cpu_count()
pool = mp.Pool(processes=cpunum-1,maxtasksperchild=1000)
pool.map(run_sim,range(sim_runs))
pool.close()
pool.join()

