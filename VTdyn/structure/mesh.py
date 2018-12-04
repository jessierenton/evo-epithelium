import numpy as np
from scipy.spatial import Delaunay, Voronoi, voronoi_plot_2d, ConvexHull
import copy

class Geometry(object):
    """Abstract Geometry object needed for Mesh."""
    
    def periodise(self,r):
        """returns coordinate r accounting for bc's"""
        raise NotImplementedError()
    
    def periodise_list(self,r):
        """returns list of coords accounting for bc's"""
        raise NotImplementedError()
         
    def retriangulate(self,centres,N_mesh):
        """Takes coordinates of set of points (centres) and number of points as arguments
        Performs Voronoi Tessellation (if calculating cell areas) or Delaunay Triangulation.
        Returns: 
            neighbours: N-dim list of (k,) arrays giving neighbour ids of each cell (k=neighbour number), 
            distances: N-dim list of (k,) arrays giving distances between each cell and its neighbours,
            sep_vector: N-dim list of (k,2) arrays giving unit vectors between each cell and its neighbours, 
            (areas: (N,) array giving area of each cell)
        """
        raise NotImplementedError()
    
    def distance(self,r0,r1):
        """returns distance between two points"""
        raise NotImplementedError()
    
    @staticmethod
    def polygon_area(points):
        n_p = len(points)
        return 0.5*sum(points[i][0]*points[(i+1)%n_p][1]-points[(i+1)%n_p][0]*points[i][1] for i in range(n_p))
        return neighbours, distances, sep_vector, areas

class Torus(Geometry):
    """Square domain with periodic boundary conditions"""
    
    def __init__(self,width,height):
        """width and height of periodicity"""
        self.width = width
        self.height = height
        
    def periodise(self,coords):
        half_width, half_height = self.width/2., self.height/2.
        for i,L in enumerate((half_width,half_height)):
            if coords[i] >= L: coords[i] -= L*2
            elif coords[i] < -L: coords[i] += L*2
        return coords
        
    def periodise_list(self,coords):
        half_width, half_height = self.width/2., self.height/2.
        for i,L in enumerate((half_width,half_height)):
            coords[np.where(coords[:,i] >= L)[0],i] -= L*2
            coords[np.where(coords[:,i] < -L)[0],i] += L*2
        return coords
    
    def retriangulate(self,centres,N_mesh):
        width,height = self.width, self.height
        centres_3x3 = np.reshape([centres+[dx, dy] for dx in [-width, 0, width] for dy in [-height, 0, height]],(9*N_mesh,2))
        vor = Voronoi(centres_3x3)
        pairs = vor.ridge_points
        neighbours = [pairs[loc[0],1-loc[1]] for loc in (np.where(pairs==k) for k in xrange(4*N_mesh,5*N_mesh))]
        sep_vectors = [centres[i]-centres_3x3[n_cell] for i,n_cell in enumerate(neighbours)]
        distances = [np.sqrt((cell_vectors*cell_vectors).sum(axis=1)) for cell_vectors in sep_vectors]
        sep_vectors = [cell_vectors/np.repeat(cell_distances[:,np.newaxis],2,axis=1) for cell_distances,cell_vectors in zip(distances,sep_vectors)]
        neighbours = [n_set%N_mesh for n_set in neighbours] 
        areas = np.abs([polygon_area(vor.vertices[polygon]) for polygon in np.array(vor.regions)[vor.point_region][4*N_mesh:5*N_mesh]])
        return neighbours, distances, sep_vectors, areas
    
    def distance(self,r0,r1):
        delta = np.abs(r0-r1)
        delta[:,0] = np.min((delta[:,0],self.width-delta[:,0]),axis=0)
        delta[:,1] = np.min((delta[:,1],self.height-delta[:,1]),axis=0)
        return np.sqrt((delta ** 2).sum(axis=1))

class TorusNoArea(Torus):
    """same as Torus geometry but does not calculate cell areas (overides retriangulate)"""           
    def retriangulate(self,centres,N_mesh):
        width,height = self.width,self.height
        centres_3x3 = np.reshape([centres+[dx, dy] for dx in [-width, 0, width] for dy in [-height, 0, height]],(9*N_mesh,2))
        vnv = Delaunay(centres_3x3).vertex_neighbor_vertices
        neighbours = [vnv[1][vnv[0][k]:vnv[0][k+1]] for k in xrange(4*N_mesh,5*N_mesh)]
        sep_vectors = [centres[i]-centres_3x3[n_cell] for i,n_cell in enumerate(neighbours)]
        distances = [np.linalg.norm(cell_vectors,axis=1) for cell_vectors in sep_vectors]
        sep_vectors = [cell_vectors/np.repeat(cell_distances[:,np.newaxis],2,axis=1) for cell_distances,cell_vectors in zip(distances,sep_vectors)]
        neighbours = [n_set%N_mesh for n_set in neighbours] 

        return neighbours,distances,sep_vectors


class Mesh(object):
    """ 
    keeps track of cell positions and neighbour relations and defines methods for moving cells and updating
    Attributes: N_cells = number of cells 
                centres = array of (x,y) values for both cell and ghost node positions
                geometry = Geometry object, e.g. Torus
                neighbours, distances, unit_vecs, areas (see Geometry class)
    """
   
    def __init__(self,centres,geometry):
        """Parameters:
        centres: (N,2) array floats
            positions of cells
        geometry: Geometry object 
        """
        self.N_mesh = len(centres)
        self.centres = centres
        self.geometry = geometry
        self.neighbours,self.distances,self.unit_vecs, self.areas = self.retriangulate()
    
    def __len__(self):
        return self.N_mesh
    
    def copy(self):
        """create a copy of Mesh object"""
        meshcopy = copy.copy(self)
        meshcopy.centres = copy.copy(meshcopy.centres)
        return meshcopy
    
    def update(self):
        """recalculate and define mesh attributes"""
        self.N_mesh = len(self.centres)
        self.neighbours, self.distances, self.unit_vecs, self.areas = self.retriangulate()
        
    def retriangulate(self):
        return self.geometry.retriangulate(self.centres,self.N_mesh)
        
    def move(self, i, dr):
        """move cell i by dr"""
        self.centres[i] = self.geometry.periodise(self.centres[i]+dr)
    
    def move_all(self, dr_array):
        """move all N cells by vectors given by (N,2) dr_array """
        self.centres = self.geometry.periodise_list(self.centres + dr_array)
    
    def add(self,pos):
        """add new cell centre"""
        self.centres = np.append(self.centres,pos,0)
    
    def remove(self,i):
        """remove cell centre"""
        self.centres = np.delete(self.centres,i,0)
        
    def voronoi(self):
        return Voronoi(self.centres)
    
    def convex_hull(self):
        return ConvexHull(self.centres())  

    def delaunay(self):
        return Delaunay(self.centres)
        
    def get_distances(self,i):
        """get distances between cell i and its neighbours"""
        return self.geometry.distance(self.centres[i],self.centres)
    
        
class MeshNoArea(Mesh):
    def __init__(self,centres,geometry):
        self.N_mesh = len(centres)
        self.centres = centres
        self.geometry = geometry
        self.neighbours,self.distances,self.unit_vecs = self.retriangulate()
    
    def update(self):
        self.N_mesh = len(self.centres)
        self.neighbours, self.distances, self.unit_vecs = self.retriangulate()

    def retriangulate(self):
        return self.geometry.retriangulate(self.centres,self.N_mesh)