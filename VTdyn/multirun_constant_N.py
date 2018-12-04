from multiprocessing import Pool #parallel processing
from itertools import repeat
import sys
from functools import partial

import os
import numpy as np
import libs.data as data
from libs.pd_lib import run_simulation,simulation_decoupled_update,prisoners_dilemma_averaged,prisoners_dilemma_accumulated
from functools import partial

"""command line arguments
1. str (av or acc): game_str. determines payoff accounting is averaged or accumulated
2, 3. ints: start_batch and end_batch. indices (run batches of 1000 simulations)
4, ... array floats: b_vals. values of b (prisoner's dilemma param) to run simulations for 
"""
game_str = sys.argv[1]
if game_str == 'acc': game = prisoners_dilemma_accumulated
elif game_str == 'av': game = prisoners_dilemma_averaged
else: raise ValueError('invalid game string')

start_batch,end_batch = int(sys.argv[2]),int(sys.argv[3])
runs_per_batch = 1000

b_vals = np.array(sys.argv[4:],dtype=float)
c,DELTA = 1.0,0.025 #prisoner dilemma params

l = 10 #population size N=lxl 
timend = 10000. #time (hours) after which simulation ends if no fixation
timestep = 12.0 #state saved every 12 hours

rand = np.random.RandomState()

outdir = 'VTpd_%s_decoupled'%(game_str)
if not os.path.exists(outdir): # if the outdir doesn't exist create it
     os.makedirs(outdir)
with open(outdir+'/info','w') as f:
    f.write('N=%d, c=%.1f, delta=%.3f'%(l*l,c,DELTA))

def run_parallel(b,i):
    """run a single simulation using simulation_decoupled_update routine indexed by i returning 1 if resulted in mutant fixation and 0 if resident fixation.
    If no fixation returns -1 and saves number of mutants at each timestep to file
    """
    rand=np.random.RandomState()
    history = run_simulation(simulation_decoupled_update,l,timestep,timend,rand,DELTA,prisoners_dilemma_averaged,(b,c),save_areas=False)
    if 0 not in history[-1].properties['type']:
        fix = 1  
    elif 1 not in history[-1].properties['type']:
        fix = 0
    else: 
        fix = -1
        data.save_N_mutant_type(history,outdir+'/incomplete_b%.1f'%b,i)
    return fix

pool = Pool(maxtasksperchild=1000) # creating a pool of workers to run simulations in parallel

"""for each b we run 1000x(end_batch-start_batch simulations). At the end of each batch of 1000 simulations write to 
file how many are fixed, lost and incomplete
"""
for b in b_vals:
    fix_results = open(outdir+'/fix%.2f'%b,'a',0)
    for i in range(start_batch,end_batch):
        text = '\r running batch %d of %d'%(i+1,end_batch)
        sys.stdout.write(text)
        sys.stdout.flush()
        fixation = np.array([f for f in pool.imap(partial(run_parallel,b),range(i*runs_per_batch,(i+1)*runs_per_batch))])  # mapping of all the calls necessary into the calling function (run parallel)
        fixed = len(np.where(fixation==1)[0])
        lost = len(np.where(fixation==0)[0])
        incomplete = len(np.where(fixation==-1)[0])
        fix_results.write('%d    %d    %d\n'%(fixed,lost,incomplete))
