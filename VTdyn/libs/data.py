import os
import numpy as np

#library of functions for saving data from a history object (list of tissues)

def save_mean_area(history,outdir,index=0):
    """saves mean area of cells in each tissue"""
    if not os.path.exists(outdir): # if the folder doesn't exist create it
         os.makedirs(outdir)
    filename = '%s/area_mean_%d'%(outdir,index)
    np.savetxt(filename,[np.mean(tissue.mesh.areas) for tissue in history])

def save_areas(history,outdir,index=0):
    """saves all areas of cells in each tissue"""
    if not os.path.exists(outdir): # if the folder doesn't exist create it
         os.makedirs(outdir)
    filename = '%s/areas_%d'%(outdir,index)
    wfile = open(filename,'w')
    for tissue in history:
        for area in tissue.mesh.areas:        
            wfile.write('%.3e    '%area)
        wfile.write('\n')
    
def save_force(history,outdir,index=0):
    """saves mean magnitude of force on cells in each tissue"""
    if not os.path.exists(outdir): # if the folder doesn't exist create it
         os.makedirs(outdir)
    wfile = open('%s/%s_%d'%(outdir,'force',index),'w')
    for tissue in history:        
        wfile.write('%.3e \n'%np.mean(np.sqrt((tissue.Force(tissue)**2).sum(axis=1))))
    wfile.close() 

def save_neighbour_distr(history,outdir,index=0,snap=-1):
    """save neighbour distributions in each tissue"""
    if not os.path.exists(outdir): # if the folder doesn't exist create it
         os.makedirs(outdir)
    wfilename = '%s/%s_%d'%(outdir,'neigh_distr',index) 
    np.savetxt(wfilename,[np.bincount([len(tissue.mesh.neighbours[i]) for i in range(len(tissue))],minlength=18) for tissue in history],fmt=(['%d']*18))
        
def save_N_cell(history,outdir,index=0):
    """save number of cells in each tissue"""
    if not os.path.exists(outdir): # if the folder doesn't exist create it
         os.makedirs(outdir)
    wfilename = '%s/%s_%d'%(outdir,'N_cell',index)  
    np.savetxt(wfilename,[len(tissue) for tissue in history],fmt=('%d'))

def save_N_mutant(history,outdir,index=0):
    """saves number of mutants in each tissue given by 'mutant' property"""
    if not os.path.exists(outdir): # if the folder doesn't exist create it
         os.makedirs(outdir)
    wfilename = '%s/%s_%d'%(outdir,'N_mutant',index)  
    np.savetxt(wfilename,[sum(tissue.properties['mutant']) for tissue in history],fmt=('%d'))

def save_N_mutant_type(history,outdir,index=0):
    """saves number of mutants in each tissue given by 'type' property"""
    if not os.path.exists(outdir): # if the folder doesn't exist create it
         os.makedirs(outdir)
    wfilename = '%s/%s_%d'%(outdir,'N_mutant',index)  
    np.savetxt(wfilename,[sum(tissue.properties['type']) for tissue in history],fmt=('%d'))

def get_cell_history(history,cell_id):
    """generate a history for a given cell id with area at each timestep, age at each timestep and 
    fate (reproduction=1 or death=0)"""
    cell_history = {'area':[],'age':[],'fate':None}
    for tissue in history:
        if cell_id not in tissue.cell_ids and len(cell_history['area']) > 0: 
            if cell_id in tissue.mother: cell_history['fate'] = 1
            else: cell_history['fate'] = 0
            break
        elif cell_id in tissue.cell_ids:
            mesh_id = np.where(tissue.cell_ids == cell_id)[0][0]
            cell_history['area'].append(tissue.mesh.areas[mesh_id])
            cell_history['age'].append(tissue.age[mesh_id])
    return cell_history

def get_cell_histories(history,start=0):
    """generate history for all cells (see above get_cell_history)"""
    return [get_cell_history(history[start:],i) for i in range(max(history[-1].cell_ids)) if len(get_cell_history(history[start:],i)['age'])>0]

def save_age_of_death(history,outdir,index=0):
    """save cell lifetimes for each cell in history"""
    if not os.path.exists(outdir): # if the folder doesn't exist create it
         os.makedirs(outdir)
    cell_histories = np.array(get_cell_histories(history))
    has_fate = np.array([h['fate'] is not None for h in cell_histories])
    cell_histories = cell_histories[has_fate]
    fates = np.array([h['fate'] for h in cell_histories],dtype=bool)
    final_age_d = [cell['age'][-1] for cell in cell_histories[fates]]
    final_age_a = [cell['age'][-1] for cell in cell_histories[~fates]]
    np.savetxt('%s/division_age_%d'%(outdir,index),final_age_d)
    np.savetxt('%s/apoptosis_age_%d'%(outdir,index),final_age_a)
    
def save_ages(history,outdir,index=0):
    """saves all cell ages for each tissue in history"""
    if not os.path.exists(outdir): # if the folder doesn't exist create it
         os.makedirs(outdir)
    filename = '%s/ages_%d'%(outdir,index)
    wfile = open(filename,'w')
    for tissue in history:
        for age in tissue.age:        
            wfile.write('%.3e    '%age)
        wfile.write('\n')    
    
def save_mean_age(history,outdir,index=0):
    """save mean age of cells for each tissue in history"""
    if not os.path.exists(outdir): # if the folder doesn't exist create it
         os.makedirs(outdir)
    filename = '%s/age_mean_%d'%(outdir,index)
    np.savetxt(filename,[np.mean(tissue.age) for tissue in history])
        
    
def save_all(history,outdir,index=0):
    save_N_cell(history,outdir,index)
    save_N_mutant(history,outdir,index)
    save_ages(history,outdir,index)
    save_neighbour_distr(history,outdir,index)
    save_force(history,outdir,index)
    try: save_areas(history,outdir,index)
    except: pass
    save_age_of_death(history,outdir,index)
    
    
    
    