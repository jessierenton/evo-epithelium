import sys
import numpy as np
import itertools
import structure
from structure.global_constants import T_D,dt
from structure.cell import Tissue, BasicSpringForceNoGrowth
import structure.initialisation as init

def print_progress(step,N_steps):
    sys.stdout.write("\r %.2f %%"%(step*100/N_steps))
    sys.stdout.flush() 
 
def run(tissue_original,simulation,N_step,skip):
    """run a given simulation for N_step iterations
    returns list of tissue objects at intervals given by skip"""
    return [tissue_original.copy()]+[tissue.copy() for tissue in itertools.islice(simulation,skip-1,N_step,skip)]

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------POISSON-CONSTANT-POP-SIZE-AND-FITNESS------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def prisoners_dilemma_averaged(cell_type,neighbour_types,b,c):
    """calculate average payoff for single cell"""
    return -c*cell_type+b*np.sum(neighbour_types)/len(neighbour_types)

def prisoners_dilemma_accumulated(cell_type,neighbour_types,b,c):
    """calculate accumulated payoff for single cell"""
    return -c*cell_type*len(neighbour_types)+b*np.sum(neighbour_types)

def get_fitness(cell_type,neighbour_types,DELTA,game,game_constants):
    """calculate fitness of single cell"""
    return 1+DELTA*game(cell_type,neighbour_types,*game_constants)

def recalculate_fitnesses(neighbours_by_cell,types,DELTA,game,game_constants):
    """calculate fitnesses of all cells"""
    return np.array([get_fitness(types[cell],types[neighbours],DELTA,game,game_constants) for cell,neighbours in enumerate(neighbours_by_cell)])

# def simulation_with_mutation_ancestor_tracking(tissue,dt,N_steps,stepsize,rand,DELTA,game,constants,initial=False):
#     """simulation loop for decoupled update rule with mutation. tracks ancestors in tissue.properties['ancestor']"""
#     mutation_rate = constants[0]
#     game_constants = constants[1:]
#     step = 0.
#     complete = False
#     while initial or not complete:
#         N= len(tissue)
#         properties = tissue.properties
#         mesh = tissue.mesh
#         step += 1
#         mesh.move_all(tissue.dr(dt))
#         if rand.rand() < (1./T_D)*N*dt:
#             fitnesses = recalculate_fitnesses(tissue.mesh.neighbours,properties['type'],DELTA,game,game_constants)
#             mother = np.where(np.random.multinomial(1,fitnesses/sum(fitnesses))==1)[0][0]
#             tissue.add_daughter_cells(mother,rand)
#             r = rand.rand()
#             if r < mutation_rate**2: properties['type'] = np.append(properties['type'],rand.randint(0,2,2))
#             elif r < mutation_rate: properties['type'] = np.append(properties['type'],[properties['type'][mother],rand.randint(2)])
#             else: properties['type'] = np.append(properties['type'],[properties['type'][mother]]*2)
#             properties['ancestor'] = np.append(properties['ancestor'],[properties['ancestor'][mother]]*2)
#             tissue.remove(mother)
#             tissue.remove(rand.randint(N)) #kill random cell
#         tissue.update(dt)
#         complete = (1 not in tissue.properties['type'] or 0 not in tissue.properties['type']) and step%stepsize==0
#         yield tissue
#
# def simulation_with_mutation(tissue,dt,N_steps,stepsize,rand,DELTA,game,constants,initial=False):
#     """simulation loop for decoupled update rule with mutation"""
#     mutation_rate = constants[0]
#     game_constants = constants[1:]
#     step = 0.
#     complete = False
#     while initial or not complete:
#         N= len(tissue)
#         properties = tissue.properties
#         mesh = tissue.mesh
#         step += 1
#         mesh.move_all(tissue.dr(dt))
#         if rand.rand() < (1./T_D)*N*dt:
#             fitnesses = recalculate_fitnesses(tissue.mesh.neighbours,properties['type'],DELTA,game,game_constants)
#             mother = np.where(np.random.multinomial(1,fitnesses/sum(fitnesses))==1)[0][0]
#             tissue.add_daughter_cells(mother,rand)
#             r = rand.rand()
#             if r < mutation_rate**2: properties['type'] = np.append(properties['type'],rand.randint(0,2,2))
#             elif r < mutation_rate: properties['type'] = np.append(properties['type'],[properties['type'][mother],rand.randint(2)])
#             else: properties['type'] = np.append(properties['type'],[properties['type'][mother]]*2)
#             tissue.remove(mother)
#             tissue.remove(rand.randint(N)) #kill random cell
#         tissue.update(dt)
#         complete = (1 not in tissue.properties['type'] or 0 not in tissue.properties['type']) and step%stepsize==0
#         yield tissue


def simulation_decoupled_update(tissue,dt,N_steps,stepsize,rand,DELTA,game,game_constants,til_fix=False):
    """simulation loop for decoupled update rule"""
    step = 0.
    complete = False
    while not til_fix or not complete:
        N= len(tissue)
        properties = tissue.properties
        mesh = tissue.mesh
        step += 1
        mesh.move_all(tissue.dr(dt))
        if rand.rand() < (1./T_D)*N*dt:
            fitnesses = recalculate_fitnesses(tissue.mesh.neighbours,properties['type'],DELTA,game,game_constants)
            mother = np.where(np.random.multinomial(1,fitnesses/sum(fitnesses))==1)[0][0]   
            tissue.add_daughter_cells(mother,rand)
            properties['type'] = np.append(properties['type'],[properties['type'][mother]]*2)
            tissue.remove(mother)
            tissue.remove(rand.randint(N)) #kill random cell
        tissue.update(dt)
        complete = (1 not in tissue.properties['type'] or 0 not in tissue.properties['type']) and step%stepsize==0  
        yield tissue

def simulation_death_birth(tissue,dt,N_steps,stepsize,rand,DELTA,game,game_constants,til_fix=False):
    """simulation loop for death-birth update rule"""
    step = 0.
    complete = False
    while not til_fix or not complete:
        N= len(tissue)
        properties = tissue.properties
        mesh = tissue.mesh
        step += 1
        mesh.move_all(tissue.dr(dt))
        if rand.rand() < (1./T_D)*N*dt:
            dead_cell = rand.randint(N)
            dead_cell_neighbours = tissue.mesh.neighbours[dead_cell]
            neighbours_by_cell = [tissue.mesh.neighbours[dcn] for dcn in dead_cell_neighbours]
            fitnesses = np.array([get_fitness(tissue.properties['type'][cell],tissue.properties['type'][neighbours],DELTA,game,game_constants) for cell,neighbours in zip(dead_cell_neighbours,neighbours_by_cell)])
            mother = rand.choice(dead_cell_neighbours,p=fitnesses/sum(fitnesses))
            tissue.add_daughter_cells(mother,rand)
            properties['type'] = np.append(properties['type'],[properties['type'][mother]]*2)
            tissue.remove(mother)
            tissue.remove(dead_cell) #kill random cell
        tissue.update(dt)
        complete = (1 not in tissue.properties['type'] or 0 not in tissue.properties['type']) and step%stepsize==0  
        yield tissue

def simulation_no_division(tissue,dt,N_steps,rand):
    """run tissue simulation with no death or division"""
    step = 0.
    while True:
        N= len(tissue)
        mesh = tissue.mesh
        step += 1
        mesh.move_all(tissue.dr(dt))
        tissue.update(dt)
        yield tissue
        
# def simulation_death_birth_radius(tissue,dt,N_steps,stepsize,rand,DELTA,game,game_constants,division_radius,initial=False):
#     """simulation loop for death-birth update rule with interactions occurring within a given radius rather than between neighbours"""
#     step = 0.
#     complete = False
#     while initial or not complete:
#         N= len(tissue)
#         properties = tissue.properties
#         mesh = tissue.mesh
#         step += 1
#         mesh.move_all(tissue.dr(dt))
#         if rand.rand() < (1./T_D)*N*dt:
#             dead_cell = rand.randint(N)
#             cell_in_range = np.where((mesh.get_distances(dead_cell)<division_radius))[0]
#             cell_in_range = cell_in_range[cell_in_range!=dead_cell]
#             neighbours_by_cell = [tissue.mesh.neighbours[cell] for cell in cell_in_range]
#             fitnesses = np.array([get_fitness(tissue.properties['type'][cell],tissue.properties['type'][neighbours],DELTA,game,game_constants)
#                                 for cell,neighbours in zip(cell_in_range,neighbours_by_cell)])
#             mother = rand.choice(cell_in_range,p=fitnesses/sum(fitnesses))
#             tissue.add_daughter_cells(mother,rand)
#             properties['type'] = np.append(properties['type'],[properties['type'][mother]]*2)
#             tissue.remove(mother)
#             tissue.remove(dead_cell) #kill random cell
#         tissue.update(dt)
#         complete = (1 not in tissue.properties['type'] or 0 not in tissue.properties['type']) and step%stepsize==0
#         yield tissue

    
def run_simulation(simulation,N,timestep,timend,rand,DELTA,game,constants,til_fix=True,save_areas=False,tissue=None):
    """initialise tissue with NxN cells and run given simulation with given game and constants.
            starts with single cooperator
            ends at time=timend OR if til_fix=True when population all cooperators (type=1) or defectors (2)
        returns history: list of tissue objects at time intervals given by timestep
            """
    if tissue is None:
        tissue = init.init_tissue_torus(N,N,0.01,BasicSpringForceNoGrowth(),rand,save_areas=False)
    tissue.properties['type'] = np.zeros(N*N,dtype=int)
    tissue.age = np.zeros(N*N,dtype=float)
    tissue = run(tissue, simulation(tissue,dt,10./dt,timestep/dt,rand,DELTA,game,constants,til_fix=False),10./dt,timestep/dt)[-1]
    tissue.properties['type'][rand.randint(N*N,size=1)]=1
    history = run(tissue, simulation(tissue,dt,timend/dt,timestep/dt,rand,DELTA,game,constants,til_fix),timend/dt,timestep/dt)
    return history
