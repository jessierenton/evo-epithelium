import numpy as np
import libs.pd_lib as lib #library for simulation routines
import libs.data as data
import libs.plot as vplt #plotting library
from structure.global_constants import *
import structure.initialisation as init
from structure.cell import Tissue, BasicSpringForceNoGrowth

"""run a single voronoi tessellation model simulation"""

l = 10 # population size N=l*l
timend = 10 # simulation time (hours)
timestep = 1. # time intervals to save simulation history

rand = np.random.RandomState()
b,c,DELTA = 10.,1.0,0.0 #prisoner's dilemma game parameters

simulation = lib.simulation_decoupled_update  #simulation routine imported from lib
game = lib.prisoners_dilemma_averaged #game imported from lib
game_parameters = (b,c)

history = lib.run_simulation(simulation,l,timestep,timend,rand,DELTA,game,game_parameters)
                 