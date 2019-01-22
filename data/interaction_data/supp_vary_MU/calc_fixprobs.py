import seaborn as sns
import numpy as np
import os

N = 100
DELTA = 0.025
c = 1.

def load_data(folder,dtype=float):    
    return [np.loadtxt(folder+'/n_%d'%n,dtype=dtype) for n in range(1,100)]

def save_mean_Lambda_CC(data_CC,str_MU):
    """saves mean Lambda_CC vals in file"""
    coop_int = np.array([np.mean(dat_n)/(i+1.) for i,dat_n in enumerate(data_CC)])
    np.savetxt('MU-%s/Lambda_CC'%str_MU,coop_int)
    return coop_int

# def fix_prob_1(b,c,Lambda_CC,DELTA):
#     """calculates fix probability for single initial cooperator"""
#     f = lambda j: -c/N +b/N*(Lambda_CC[j-1]*N-j)/(N-j)
#     return 1./N +DELTA/N*sum(sum(f(j) for j in range(1,k+1)) for k in range(1,N))

def fix_prob_1(b,c,Lambda_CC,DELTA):
    """calculates fix probability for single initial cooperator"""
    f = lambda j: (Lambda_CC[j-1]-float(j)/N)/float(N-j)
    return 1./N +DELTA/N*( (-c*(N-1.))/2. + b*sum(sum(f(j) for j in range(1,k+1)) for k in range(1,N)))


def critical_ratio(str_MU):
    try: Lambda_CC = np.loadtxt('MU-%s/Lambda_CC'%str_MU)
    except: 
        data_CC = load_data('MU-%s/wints_CC'%str_MU)
        Lambda_CC = save_mean_Lambda_CC(data_CC,str_MU)
    f = lambda j: (Lambda_CC[j-1]-float(j)/N)/float(N-j)
    return (N-1.)/2. *1/(sum(sum(f(j) for j in range(1,k+1)) for k in range(1,N)))

def save_critical_ratios(str_MU_list,fname):
    with open(fname,'w') as f:
        for str_MU in str_MU_list:
            f.write('%s    %.3f \n'%(str_MU,critical_ratio(str_MU)))

def save_fix_probs(DELTA_vals,b_vals,fname,str_MU):
    try: Lambda_CC = np.loadtxt('MU-%s/Lambda_CC'%str_MU)
    except: 
        data_CC = load_data('MU-%s/wints_CC'%str_MU)
        Lambda_CC = save_mean_Lambda_CC(data_CC,str_MU)
    with open(fname,'w') as f:
        for b in b_vals:
            f.write('%.1f     %.3f\n'%(b,fix_prob_1(b,c,Lambda_CC,DELTA)*100))

str_MU_list = ('0.1','1','10','25','50','100','250')
for str_MU in str_MU_list:
    b_vals = np.arange(1.5,4.5,0.5)
    save_fix_probs(DELTA,b_vals,'MU-%s_fixprobs'%str_MU,str_MU)
    
save_critical_ratios(str_MU_list,'critical_ratios')