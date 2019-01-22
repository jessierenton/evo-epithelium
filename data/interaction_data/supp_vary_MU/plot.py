import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
import seaborn as sns
import os

def read_data(filename):
    dat = np.loadtxt(filename,dtype=float).T
    fix = dat[0].sum()
    lost = dat[1].sum()
    return fix,lost
    
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
    if legend: plt.legend(loc='best',fontsize=10,frameon=False)

def plot_neutral_line(ax,start,stop):
    plt.plot((start,stop),[1.0]*2,ls=':',color='grey',label=None)

def plot_theory(filename,ax,color,start_index=None,end_index=None,label=None,percent=False,zorder=None):
    fix_theory = np.loadtxt(filename).T
    if not percent: ax.plot(fix_theory[0][start_index:end_index],fix_theory[1][start_index:end_index]*100,color=color,label=label,zorder=zorder)
    else: ax.plot(fix_theory[0][start_index:end_index],fix_theory[1][start_index:end_index],color=color,label=label,zorder=zorder)

def plot_critical_ratios(filename,ax):
    data = np.loadtxt(filename).T
    ax.plot(data[0],data[1],ls='-',marker='o')

sns.set_style('white')    

cp = sns.color_palette('colorblind')
sns.set_palette(cp)

#Supplementary information Figure 2
fig,ax = plt.subplots()
plot_neutral_line(ax,2,4)
plot_theory('MU-0.1_fixprobs',ax,cp[0],percent=True,start_index=1,label=r'$\mu=0.1$')
plot_theory('MU-1_fixprobs',ax,cp[1],percent=True,start_index=1,label=r'$\mu=1$')
plot_theory('MU-10_fixprobs',ax,cp[2],percent=True,start_index=1,label=r'$\mu=10$')
plot_theory('MU-25_fixprobs',ax,cp[3],percent=True,start_index=1,label=r'$\mu=25$')
plot_theory('MU-50_fixprobs',ax,cp[4],percent=True,start_index=1,label=r'$\mu=50$')
plot_theory('MU-100_fixprobs',ax,cp[5],percent=True,start_index=1,label=r'$\mu=100$')
plot_theory('MU-250_fixprobs',ax,cp[6],percent=True,start_index=1,label=r'$\mu=250$')
formatting(r'Benefit-to-cost ratio, $b/c$',r'Fixation Probability, $\rho _C$ (%)',large=False,legend=True)
plt.xticks(np.arange(2,4.5,0.5))
plt.savefig('vary_MU_fixprob.pdf')

plt.show()

#Supplementary information Figure 3
fig,ax = plt.subplots()
plot_critical_ratios('critical_ratios',ax)
ax.annotate(r'$(b/c)^*=2.78$',xy=(50,2.77),xytext=(45,2.65),arrowprops=dict(facecolor='black',arrowstyle='simple'))
formatting(r'$\mu$',r'$(b/c)^*$',large=False,legend=True)
plt.savefig('critical_ratios.pdf')
plt.show()

