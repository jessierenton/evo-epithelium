import numpy as np
import copy
from functools import partial
import global_constants as gc
from global_constants import EPS, L0, MU, ETA
              
class Tissue(object):    
    
    """Defines a tissue comprised of cells which can move, divide and be extruded"""
    
    def __init__(self,mesh,force,cell_ids,next_id,age,mother,properties=None):
        """ Parameters:
        mesh: Mesh object
            defines cell locations and neighbour connections
        force: Force object
            defines force law between neighbouring cells
        cell_ids: (N,) array ints
             unique id for each cell (N is number of cells)
        next_id: int
            next available cell id
        age: (N,) array floats
            age of each cell
        mother: (N,) array ints
            id of mother for each cell (-1 for initial cells)
        properties: dict or None
            dictionary available for any other cell properties
            
        """
        self.mesh = mesh
        self.Force = force
        N = len(mesh)
        self.cell_ids = cell_ids
        self.next_id = next_id
        self.age = age
        self.mother = mother
        self.properties = properties or {}
        
    def __len__(self):
        return len(self.mesh)
    
    def copy(self):
        """create a copy of Tissue"""
        return Tissue(self.mesh.copy(),self.Force,self.cell_ids.copy(),self.next_id,self.age.copy(),self.mother.copy(),self.properties.copy())
            
    def update(self,dt):
        self.mesh.update()
        self.age += dt
    
    def remove(self,idx_list):
        """remove a cell from tissue"""
        self.mesh.remove(idx_list)
        self.cell_ids = np.delete(self.cell_ids,idx_list)
        self.age = np.delete(self.age,idx_list)
        self.mother = np.delete(self.mother,idx_list)
        for key,val in self.properties.iteritems():
            self.properties[key] = np.delete(val,idx_list)
        
    def add_daughter_cells(self,i,rand):
        """add pair of new cells after a cell division"""
        angle = rand.rand()*np.pi
        dr = np.array((EPS*np.cos(angle),EPS*np.sin(angle)))
        new_cen1 = self.mesh.centres[i] + dr
        new_cen2 = self.mesh.centres[i] - dr
        self.mesh.add([new_cen1,new_cen2])
        self.cell_ids = np.append(self.cell_ids,[self.next_id,self.next_id+1])
        self.age = np.append(self.age,[0.0,0.0])
        self.mother = np.append(self.mother,[self.cell_ids[i]]*2)
        self.next_id += 2
        
    def dr(self,dt): 
        """calculate distance cells move due to force law in time dt"""  
        return (dt/ETA)*self.Force(self)
        
        
class Force(object):
    """Abstract force object"""
    def force(self):
        """returns (N,2) array floats giving vector force on each cell"""
        raise NotImplementedError()
    
    def magnitude(self,tissue):
        """returns (N,) array floats giving magnitude of force on each cell"""
        return np.sqrt(np.sum(np.sum(self.force(tissue),axis=0)**2))
    
    def __call__(self, tissue):
        return self.force(tissue)
        
class BasicSpringForceTemp(Force):
    
    def __init__(self,MU=MU):
        self.MU=MU
    
    def force(self,tissue):
        return np.array([self.force_i(tissue,i,dist,vec,neigh) for i,(dist,vec,neigh) in enumerate(zip(tissue.mesh.distances,tissue.mesh.unit_vecs,tissue.mesh.neighbours))])
    
    def force_i(self):
        """returns force on cell i"""
        raise Exception('force law undefined')

class BasicSpringForceNoGrowth(BasicSpringForceTemp):
    
    def force_i(self,tissue,i,distances,vecs,n_list):
        if tissue.age[i] >= 1.0 or tissue.mother[i] == -1: pref_sep = L0
        else: pref_sep = (tissue.mother[n_list]==tissue.mother[i])*((L0-EPS)*tissue.age[i]+EPS-L0) +L0
        return (self.MU*vecs*np.repeat((distances-pref_sep)[:,np.newaxis],2,axis=1)).sum(axis=0)    


