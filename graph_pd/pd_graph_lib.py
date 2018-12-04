# from igraph import Graph
import numpy as np
import itertools
import multiprocessing as mp

#Library of functions for running EGT simulations with a prisoners dilemma game determining payoffs and a Moran process with either birth-death or death-birth update rules

def prisoners_dilemma_averaged(cell_type,neighbour_types,b,c):
    """returns payoff to cell calculated according to pd with averaged accounting"""
    return -c*cell_type+b*np.sum(neighbour_types)/len(neighbour_types)

def prisoners_dilemma_accumulated(cell_type,neighbour_types,b,c):
    """returns payoff to cell calculated according to pd with accumulated accounting"""
    return np.sum(cell_type*neighbour_types*(b-c)+cell_type*(1-neighbour_types)*-c+(1-cell_type)*neighbour_types*b)
    
def get_all_fitnesses(cell_types,neighbours,DELTA,game,game_params):
    """returns list of cell fitnesses"""
    return 1+DELTA*np.array([game(cell_types[cell],cell_types[cell_neighbours],*game_params) for cell,cell_neighbours in enumerate(neighbours)])

def get_fitness(cell_type,cell_neighbour_types,DELTA,game,game_params):
    """returns fitness of single cell"""
    return 1+DELTA*game(cell_type,cell_neighbour_types,*game_params)
    
def death_birth_update(N,cell_types,neighbours,DELTA,game,game_params,rand,til_fix=False):
    """implement death-birth update and yield new list cell types"""
    while not til_fix or (0 in cell_types and 1 in cell_types):
        dead_cell = rand.randint(N)
        dead_cell_neighbours = neighbours[dead_cell]
        neighbour_fitnesses = [get_fitness(cell_types[cell],cell_types[cell_neighbours],DELTA,game,game_params) for cell,cell_neighbours in zip(dead_cell_neighbours,[neighbours[dc_n] for dc_n in dead_cell_neighbours])]
        mother_cell = rand.choice(dead_cell_neighbours,p=neighbour_fitnesses/np.sum(neighbour_fitnesses))
        cell_types[dead_cell] = cell_types[mother_cell]
        yield cell_types

def birth_death_update(N,cell_types,neighbours,DELTA,game,game_params,rand,til_fix=False):
    """implement birth-death update and yield new list cell types"""
    while not til_fix or (0 in cell_types and 1 in cell_types):
        fitnesses = get_all_fitnesses(cell_types,neighbours,DELTA,game,game_params)
        mother_cell = rand.choice(N,p=fitnesses/np.sum(fitnesses))
        dead_cell = rand.choice(neighbours[mother_cell])
        cell_types[dead_cell] = cell_types[mother_cell]
        yield cell_types

def death_birth_update_neutral_clones(N,cell_types,neighbours,rand):
    """implement death-birth update when all cells have equal fitness"""
    while (np.max(cell_types)!=np.min(cell_types)):
        dead_cell = rand.randint(N)
        dead_cell_neighbours = neighbours[dead_cell]
        mother_cell = rand.choice(dead_cell_neighbours)
        cell_types[dead_cell] = cell_types[mother_cell]
        yield cell_types

def death_birth_update_migration(N,cell_types,neighbours,DELTA,game,game_params,migration_strength,rand,til_fix=False):
    """implement death-birth update with migration"""
    while not til_fix or (0 in cell_types and 1 in cell_types):
        # death-birth event
        dead_cell = rand.randint(N)
        dead_cell_neighbours = neighbours[dead_cell]
        neighbour_fitnesses = [get_fitness(cell_types[cell],cell_types[cell_neighbours],DELTA,game,game_params) for cell,cell_neighbours in zip(dead_cell_neighbours,[neighbours[dc_n] for dc_n in dead_cell_neighbours])]
        mother_cell = rand.choice(dead_cell_neighbours,p=neighbour_fitnesses/np.sum(neighbour_fitnesses))
        cell_types[dead_cell] = cell_types[mother_cell]
        #migration as a swap
        if rand.rand() < migration_strength:
            swap1 = rand.randint(N)
            swap2 = rand.choice(neighbours[swap1])
            cell_types[swap1],cell_types[swap2] = cell_types[swap2],cell_types[swap1]        
        yield cell_types

def run(cell_types_t0,update,N_step,skip):
    return [cell_types_t0]+[cell_types.copy() for cell_types in itertools.islice(update,skip-1,N_step,skip)]

def run_death_birth_simulation_til_fixation(neighbours,DELTA,game,game_params,timend,timestep,rand,return_fix = False,migration_strength=None):
    N = len(neighbours)
    cell_types = np.zeros(N)
    cell_types[rand.randint(N)]=1
    if migration_strength is None: history = run(cell_types.copy(),death_birth_update(N,cell_types,neighbours,DELTA,game,game_params,rand,til_fix=True),timend,timestep)
    else: history = run(cell_types.copy(),death_birth_update_migration(N,cell_types,neighbours,DELTA,game,game_params,migration_strength,rand,til_fix=True),timend,timestep)
    if return_fix:
        number_mutants = sum(history[-1])
        if number_mutants >= N: fix = 1
        elif number_mutants <= 0: fix = 0
        else: fix =  -1
        return fix, history
    else: return history
    
def run_birth_death_simulation_til_fixation(neighbours,DELTA,game,game_params,timend,timestep,rand,return_fix = False):
    N = len(neighbours)
    cell_types = np.zeros(N)
    cell_types[rand.randint(N)]=1
    history = run(cell_types.copy(),birth_death_update(N,cell_types,neighbours,DELTA,game,game_params,rand,til_fix=True),timend,timestep)
    if return_fix:
        number_mutants = sum(history[-1])
        if number_mutants >= N: fix = 1
        elif number_mutants <= 0: fix = 0
        else: fix =  -1
        return fix, history
    else: return history
    
def run_death_birth_simulation_neutral_clones(neighbours,timend,timestep,rand):
    N = len(neighbours)
    cell_types=np.arange(N)
    history = run(cell_types.copy(),death_birth_update_neutral_clones(N,cell_types,neighbours,rand),timend,timestep)
    return history

def run_death_birth_simulation(neighbours,DELTA,game,game_params,timend,timestep,rand):
    N = len(neighbours)
    cell_types = np.zeros(N)
    cell_types[rand.randint(N)]=1
    history = run(cell_types.copy(),death_birth_update(N,cell_types,neighbours,DELTA,game,game_params,rand),timend,timestep)
    return history