import seaborn as sns
import numpy as np
import os

N = 100
DELTA = 0.025
c = 1.

def save_mean_Lambda_CC(data_CC):
    """saves mean Lambda_CC vals in file"""
    coop_int = np.array([np.mean(dat_n)/(i+1.) for i,dat_n in enumerate(data_CC)])
    np.savetxt('Lambda_CC',coop_int)
    return coop_int

def fix_prob_1(b,c,Lambda_CC,DELTA):
    """calculates fix probability for single initial cooperator"""
    f = lambda j: -c/N +b/N*(Lambda_CC[j-1]*N-j)/(N-j)
    return 1./N +DELTA/N*sum(sum(f(j) for j in range(1,k+1)) for k in range(1,N))

def save_fix_probs(DELTA_vals,b_vals,fname):
    try: Lambda_CC = np.loadtxt('Lambda_CC')
    except: 
        data_CC = load_data('wints_CC')
        Lambda_CC = save_mean_Lambda_CC(data_CC)
    with open(fname,'w') as f:
        for b in b_vals:
            f.write('%.1f     %.3f\n'%(b,fix_prob_1(b,c,Lambda_CC,DELTA)*100))

if __name__ == '__main__':   
    b_vals = np.arange(1.5,8.5,0.5)
    save_fix_probs(DELTA,b_vals,'VTpd_av_decoupled_theory')