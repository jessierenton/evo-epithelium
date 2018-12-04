import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from igraph import Graph
import os
current_palette = sns.set_palette(sns.color_palette("hls", 7))

large_bool = False

def confint(p,n):
    z = 1.96
    return z*np.sqrt(p*(1-p)/n)

def formatting(xlabel,ylabel,large=False,legend=False):
    if large: labelsize,ticksize = 26,18
    else: labelsize,ticksize = 16,12   
    plt.xlabel(xlabel,size=labelsize,labelpad=10)
    plt.ylabel(ylabel,size=labelsize,labelpad=10)
    plt.xticks(size=ticksize)
    plt.yticks(size=ticksize)
    plt.tight_layout()
    if legend: plt.legend(loc='best',fontsize=10)

fig = plt.figure()
sns.set_style('white')
current_palette = sns.color_palette()

delta=0.025
b_vals = np.arange(2,12)
file_prefix = 'm'
dir_ = 'EGTpd_av_db/migration/hex_graph/'
prefix_length = len(file_prefix)
files = {f[prefix_length:prefix_length+3]:dir_+f for f in os.listdir(dir_) if f[:prefix_length]==file_prefix}
m_vals = [0.0,0.1,0.2,0.5,1.0]
plt.plot(b_vals,[1.]*len(b_vals),ls=':',color='grey')
for m in m_vals: 
    fix_data = np.loadtxt(files['%.1f'%m],skiprows=1,dtype=float)
    fix_data[:,1]= fix_data[:,1]/(fix_data[:,2]+fix_data[:,1])    
    fix_data = fix_data[:,:2].T
    plt.plot(fix_data[0],fix_data[1]*100,linestyle='',marker='o',label=r'$m=%s$'%m)
# plt.plot([6.7]*7,np.arange(7),ls='--')
plt.yticks((0.5,1.0,1.5,2.0))
legend = plt.legend(loc='upper left',fontsize=10,frameon=False)
sns.set_style('white')
sns.set_context('talk')
formatting(r'Benefit-to-cost ratio, $b/c$',r'Fixation probability, $\rho_C$ (%)',large=False)
plt.savefig('fixprobs.pdf')
plt.show()
