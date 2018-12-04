import os
import sys
import numpy as np
import itertools
import structure
from structure.global_constants import T_D,dt
from structure.cell import Tissue, BasicSpringForceNoGrowth
import structure.initialisation as init

        
def run(tissue_original,simulation,N_step,skip):
    """run a given simulation for N_step iterations
    returns list of tissue objects at intervals given by skip"""
    return [tissue_original.copy()]+[tissue.copy() for tissue in itertools.islice(simulation,skip-1,N_step,skip)]

def run_generator(simulation,N_step,skip):
    """generator for running a given simulation for N_step iterations
    returns generator for of tissue objects at intervals given by skip"""
    return itertools.islice(simulation,0,N_step,skip)

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------ SIMULATION ROUTINES ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def simulation_ancestor_tracking(tissue,dt,N_steps,stepsize,rand,til_fix=False):
    """simulation loop for neutral process tracking ancestor ids"""
    complete=False
    step = 0.
    while not til_fix or not complete:
        N= len(tissue)
        properties = tissue.properties
        mesh = tissue.mesh
        mesh.move_all(tissue.dr(dt))
        if rand.rand() < (1./T_D)*N*dt:
            mother = rand.randint(N)
            tissue.add_daughter_cells(mother,rand)
            properties['ancestor'] = np.append(properties['ancestor'],[properties['ancestor'][mother]]*2)
            tissue.remove(mother)
            tissue.remove(rand.randint(N)) #kill random cell
        tissue.update(dt)
        complete = (np.all(tissue.properties['ancestor']==tissue.properties['ancestor'][0]) and step%stepsize==0)
        step += 1 
        yield tissue

def initialise_tissue_ancestors(N,dt,timend,timestep,rand):  
    """initialise tissue and run simulation until timend returning final state"""              
    tissue = init.init_tissue_torus(N,N,0.01,BasicSpringForceNoGrowth(),rand,save_areas=False)
    tissue.properties['ancestor'] = np.arange(N*N)
    tissue.age = np.zeros(N*N,dtype=float)
    tissue = run(tissue,simulation_ancestor_tracking(tissue,dt,timend/dt,timestep/dt,rand),timend/dt,timestep/dt)[-1]
    tissue.properties['ancestor']=np.arange(N*N)
    return tissue

def run_simulation(simulation,N,timestep,timend,rand,til_fix=True,save_areas=False,tissue=None):
    """initialise tissue with NxN cells and run given simulation with given game and constants.
            starts with single cooperator
            ends at time=timend OR if til_fix=True when population all cooperators (type=1) or defectors (2)
        returns history: list of tissue objects at time intervals given by timestep
            """
    if tissue is None:
        tissue = init.init_tissue_torus(N,N,0.01,BasicSpringForceNoGrowth(),rand,save_areas=False)
        tissue.properties['type'] = np.zeros(N*N,dtype=int)
        tissue.age = np.zeros(N*N,dtype=float)
        tissue = run(tissue, simulation(tissue,dt,10./dt,timestep/dt,rand,True),10./dt,timestep/dt)[-1]
        tissue.properties['type'][rand.randint(N*N,size=1)]=1
    history = run(tissue, simulation(tissue,dt,timend/dt,timestep/dt,rand,til_fix),timend/dt,timestep/dt)
    return history

