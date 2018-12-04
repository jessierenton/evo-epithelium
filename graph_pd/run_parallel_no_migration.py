from igraph import Graph
import numpy as np
import multiprocessing as mp
from functools import partial
import os

from pd_graph_lib import run_death_birth_simulation_til_fixation, prisoners_dilemma_averaged, prisoners_dilemma_accumulated

#calculates fixation probabilities from a large number of simulation using the run_death_birth_simulation routine for a given
#graph structure  

nlist_file = 'vt_graph_6' #import neighbour list for each node in graph. possible files vt_graph_X or hex_graph
c,DELTA = 1.0,0.025 
b_vals = np.hstack((np.arange(2,4),np.arange(4,9,0.5),np.arange(9,13))) #values of b to calc fixprob for
runs=int(1e5) #number of simulations from which to calc fixprob

def read_neighbour_list(fname):
    """read in neighbours from file"""
    neighbours = []
    with open(fname,'r') as f:
        for line in f:
            neighbours.append(np.fromstring(line,dtype=int,sep=' '))
    return neighbours

def run_parallel(b,i):
    """run a single simulation with death-birth update until fixation returning +1 (mutants fixed), 0 (mutants died out), or -1 (incomplete).
        if incomplete saves history"""
    if i%int(1e3)==0: print 'b = %.1f; %d runs complete' %(b,i)
    rand = np.random.RandomState()
    fix,history = run_death_birth_simulation_til_fixation(neighbours,DELTA,game,(b,c),timend,timestep,rand,return_fix=True)
    if fix == -1: np.savetxt('%s/incomplete_b%d_%d'%(outdir,b,i),[sum(types) for types in history])
    return fix

fname = nlist_file
game = prisoners_dilemma_averaged
outdir = 'EGTpd_av_db/no_migration'
outfile = '%s/%s'%(outdir,fname)
if not os.path.exists(outdir): # if the outdir doesn't exist create it
     os.makedirs(outdir)  
neighbours = read_neighbour_list(nlist_file)
neighbours = [np.array(n,dtype=int) for n in neighbours]

N = len(neighbours)

timend,timestep = 100000.,1.0

cpunum=mp.cpu_count()
pool = mp.Pool(processes=cpunum,maxtasksperchild=1000)

fix_results = open(outfile,'w',0)
fix_results.write('#b    fixed    lost    \n')

for b in b_vals:
    fixation = np.array([f for f in pool.imap(partial(run_parallel,b),range(int(runs)))])  
    fixed = len(np.where(fixation==1)[0])
    nofixed = len(np.where(fixation==0)[0])
    incomplete = len(np.where(fixation==-1)[0])
    fix_results.write('%.1f    %d    %d    %d\n'%(b,fixed,nofixed,incomplete))

