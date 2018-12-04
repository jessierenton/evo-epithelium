#run file to produce figure 4 (Lambda_CC distribution)
#also contains functions for plotting distributions of number of (weighted) interactions
# also generates Lambda_CC file which saves mean Lambda_CC for each n

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sns.set_style('white')

def load_data(folder,dtype=float):    
    return [np.loadtxt(folder+'/n_%d'%n,dtype=dtype) for n in range(1,100)]

def save_mean_Lambda_CC(data_CC):
    """saves mean Lambda_CC vals in file"""
    coop_int = np.array([np.mean(dat_n)/(i+1.) for i,dat_n in enumerate(data_CC)])
    np.savetxt('Lambda_CC',coop_int)

def plot_mean_std(data,ax,label=None,show_error=True):
    """for a given data set plot the mean and standard deviation (if show_error=True)"""
    x = np.arange(1,100)
    mean = np.array([np.mean(data_n) for data_n in data])
    if show_error: std = np.array([np.std(data_n) for data_n in data])
    ax.plot(x,mean,label=label)
    if show_error: ax.fill_between(x,mean-std,mean+std,alpha=0.3) 
    
def plot_Lambda_CC(data,ax,label=None,show_error=True,discrete=False):
    """plot lambda_CC mean and standard deviation (if show_error=True)"""
    x=np.arange(0.01,1.,0.01)
    Lambda_CC = np.array([np.mean(data_n)/(i+1.) for i,data_n in enumerate(data)])
    if discrete and not show_error: ax.plot(x,Lambda_CC,label=label,ls='',marker='.') 
    if not discrete: ax.plot(x,Lambda_CC,label=label)
    if show_error:
        std = np.array([np.std(data_n)/(i+1.) for i,data_n in enumerate(data)])
        if discrete:
            ax.errorbar(x,Lambda_CC,yerr=std/2,label=label,ls='',marker='.',elinewidth=0.7)
        else: 
            ax.fill_between(x,Lambda_CC-std/2,Lambda_CC+std/2,alpha=0.3)  
        
def interactions_plot():
    """plots number CC and CD interactions for each cluster size n"""
    data = load_data('ints_CC'),load_data('ints_CD')
    fig,ax = plt.subplots()
    plot_mean_std(data_CC,ax,'C-C interactions')
    plot_mean_std(data_CD,ax,'C-D interactions')
    plt.xlabel('cluster size, n')
    plt.legend(loc='best')
    plt.savefig('interactions.pdf')

def weighted_interactions_plot():
    """plots number CC and CD interactions weighted by neighbour number for each cluster size n"""
    data_CC,data_CD = load_data('wints_CC'),load_data('wints_CD')
    fig,ax = plt.subplots()
    plot_mean_std(data_CC,ax,'C-C interactions (degree-weighted)')
    plot_mean_std(data_CD,ax,'C-D interactions (degree-weighted)')
    plt.xlabel('cluster size, n')
    plt.legend(loc='best')
    plt.savefig('interactions_weighted.pdf')
    
def Lambda_CC_plot(data_CC,show_error=True,discrete=False):
    """creates and saves lambda_CC plot"""
    sns.set_context('paper')
    fig,ax = plt.subplots()
    plot_Lambda_CC(data_CC,ax,show_error=show_error,discrete=discrete)
    formatting(r'Fraction of cooperators, $n/N$',r'$\Lambda^{CC}_n$',large=False)
    plt.savefig('Lambda_CC.pdf')

def formatting(xlabel,ylabel,large=False):
    if large: labelsize,ticksize = 26,18
    else: labelsize,ticksize = 16,12   
    plt.xlabel(xlabel,size=labelsize,labelpad=10)
    plt.ylabel(ylabel,size=labelsize,labelpad=10)
    plt.xticks(size=ticksize)
    plt.yticks(size=ticksize)
    plt.tight_layout() 
    
#Figure 4
data_CC = load_data('wints_CC')
save_mean_Lambda_CC(data_CC)
Lambda_CC_plot(data_CC,show_error=True,discrete=True)
    
