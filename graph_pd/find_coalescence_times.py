import numpy as np
import os

#calculates matrix of coalescence times tau_ij which give the expected meeting time of two random walks starting at vertices i and j on
#a graph. From these calculates the fixation probabilities and critical benefit-to-cost ratios for a single cooperator mutant on the graph
#plating an additive prisoners dilemma (donation) game (as derived in Allen et al 2017 "Evolutionary dynamics on any population structure")   


outdir = 'EGTpd_av_db/allen_result/'
if not os.path.exists(outdir): 
     os.makedirs(outdir)  
     
     
def load_neighbours(fname):
    """read in neighbours of each vertex from a file"""
    neighbours = []
    with open(fname,'r') as f:
        for line in f:
            neighbours.append(np.fromstring(line,dtype=int,sep=' '))
    return neighbours

def get_prob_mat(adj_matrix):
    """return 1-step probability matrix for given adjacency matrix"""
    prob_mat = adj_matrix/np.sum(adj_matrix,axis=1,keepdims=True)
    if True in np.isnan(prob_mat):
        raise Exception('Graph not connected')
    return prob_mat
    
def get_stationary_distr(adj_matrix):
    """returns stationary distribution (pi_i) for random walks on graph defined by adjacency matrix
    pi_i = k_i/K where k_i is the degree of vertex i and K=sum_i(k_i)"""
    return np.sum(adj_matrix,axis=1)/np.sum(adj_matrix)
    
def get_adjacency_from_neighbours(neighbours):
    """get adjacency matrix from a list of neighbours by vertex"""
    adj_matrix = np.zeros((len(neighbours),len(neighbours)),dtype=float)
    for vertex,vertex_neighbours in enumerate(neighbours):
        adj_matrix[vertex][vertex_neighbours] = 1.
    return adj_matrix

def build_M1(p_matrix):
    """M1 consists of nxn block matrices M1_ij s.t. [M1_ij]_lm = (p_lm-2d_lm)d_ij where d_ij is kronecker-delta"""
    n = len(p_matrix)
    M1_block = p_matrix - np.identity(n)*2
    return np.block([[M1_block if i == j else np.zeros((n,n)) for j in range(n)] for i in range(n)])
    
def build_M2(p_matrix):
    """M2 consists of nxn block matrices M2_ij s.t. [M2_ij]_lm = (p_im)(d_lj) where d_ij is kronecker-delta"""
    n = len(p_matrix)
    block = lambda i,j: np.stack([p_matrix[i] if l==j else np.zeros(n) for l in range(n)])
    return np.block([[block(i,j) for j in range(n)] for i in range(n)])

def build_M(p_matrix):
    """M is an (n^2)x(n^2) matrix. M=M1+M2 except for i(n+1)th rows for which the i(n+1)th component is 1 and all others are 0"""
    n=len(p_matrix)
    M = build_M1(p_matrix)+build_M2(p_matrix)
    for i in range(n):
        M[i*n+i,:]=0
        M[i*n+i,i*n+i]=1
    return M

def find_tau_matrix(p_matrix):
    """find coalescence times tau_ij by solving equation b = Mt
    b is an n^2-dimensional vector s.t. b_(in+j)= {-2 i!=j 
                                                  {0  i==j
    tau_in+j is an n^2-dim vector s.t. tau_in+j = tau_ij
    return tau in (nxn)-dim matrix form 
    """
    n = len(p_matrix)
    b = np.array([-2 if i!=j else 0 for i in range(n) for j in range(n)])
    M = build_M(p_matrix)
    tau = np.reshape(np.linalg.solve(M,b),(n,n))
    return tau

def find_tau_ii_plus(tau_mat,p_matrix):
    """returns tau_ii_plus calculated from tau matrix and probability matrix
    tau_ii_plus is expected remeeting time for two random walks starting at the same vertex i"""
    return np.array([1 + np.dot(p_vec,tau_vec) for p_vec,tau_vec in zip(tau_mat,p_matrix)])
    
def get_tau_1(pi_vec,tau_plus_vec):
    """tau_n is expected remeeting time for two ends of a random walk of length n started from the 
    stationary distribution. Returns tau_1"""
    return np.dot(pi_vec,tau_plus_vec) - 1.

def get_tau_2(pi_vec,tau_plus_vec):
    """as above"""
    return np.dot(pi_vec,tau_plus_vec) - 2.

def get_tau_3(pi_vec,tau_plus_vec,p_matrix):
    """as above"""
    p_ii_2_vec = np.array([np.dot(p,pT) for p,pT in zip(p_matrix,p_matrix.T)])
    return np.sum(pi_vec*tau_plus_vec*(1+p_ii_2_vec)) - 3.
    
def find_taus(tau_matrix,adj_matrix,p_matrix):
    """calculates tau_1,tau_2,tau_3 (as above)"""
    tau_plus_vec = find_tau_ii_plus(tau_matrix,p_matrix)
    pi_vec = get_stationary_distr(adj_matrix)
    
    tau_1 = get_tau_1(pi_vec,tau_plus_vec)
    tau_2 = get_tau_2(pi_vec,tau_plus_vec)
    tau_3 = get_tau_3(pi_vec,tau_plus_vec,p_matrix)
    
    return tau_1,tau_2,tau_3

def find_critical_ratio_av(tau_1,tau_2,tau_3):
    """returns critical ratio defined as the benefit to cost ratio (b/c) at which cooperation is neutral"""
    return tau_2/(tau_3-tau_1)
    
def get_fix_probs_av(b_vals,N,DELTA,c,tau_1,tau_2,tau_3):    
    """returns fixation probabilities for given set of b values"""
    return np.vstack((b_vals,1./N + DELTA/(2*N)*(-c*tau_2+b_vals*(tau_3-tau_1)))).T
    
def save_fix_probs(neighbour_file,b_vals,DELTA):   
    """prints critical benefit-to-cost ratio and saves fixprobs for given file containing adjacency matrix info"""
    adj_matrix = get_adjacency_from_neighbours(load_neighbours(neighbour_file))
    p_matrix = get_prob_mat(adj_matrix)
    tau_matrix = find_tau_matrix(p_matrix)
    taus = find_taus(tau_matrix,adj_matrix,p_matrix)
    print find_critical_ratio_av(*taus) 
    fix_probs = get_fix_probs_av(b_vals,len(p_matrix),DELTA,1.0,*taus)
    np.savetxt(outdir+neighbour_file,fix_probs,fmt=('%2.2f','%.5f'),header='b  fixprob')
    
files = ['vt_graph_%d'%i for i in range(1,7)]+['hex_graph']
b_vals = np.arange(2,13)
DELTA = 0.025

for f in files:
    save_fix_probs(f,b_vals,DELTA)
    
    
