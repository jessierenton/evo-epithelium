import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sns.set_style('white')
sns.set_palette('colorblind')

def formatting(xlabel,ylabel,large=False,legend=False):
    if large: labelsize,ticksize = 26,18
    else: labelsize,ticksize = 16,12   
    plt.xlabel(xlabel,size=labelsize,labelpad=10)
    plt.ylabel(ylabel,size=labelsize,labelpad=10)
    plt.xticks(size=ticksize)
    plt.yticks(size=ticksize)
    plt.tight_layout()
    if legend: plt.legend(loc='best',fontsize=10,frameon=False)

def load_data(data_type,MU_vals,start_index=1):
    return {MU:np.array([np.loadtxt('MU%.1f/%s_%d'%(MU,data_type,i))[start_index:] for i in range(3)]) for MU in MU_vals}
    
MU_vals = -np.array((0.1,1,10,25,50,100,250))
timestep = 1.0

def plot_force(MU_vals):
    force_data = load_data('force',MU_vals)
    t_vals = [timestep*i for i in range(len(force_data[MU_vals[0]][0]))]
    fig = plt.figure()
    for MU in MU_vals:
        plt.plot(t_vals[::10],np.mean(force_data[MU],axis=0)[::10],label=r'\mu = %d'%MU)
    plt.legend(loc='best')
    plt.show()
    
def plot_neigh_distr(MU_vals):
    neigh_data = load_data('neigh_distr',MU_vals,int(10/timestep))
    neigh_data = {MU:np.mean(data.reshape(data.shape[0]*data.shape[1],data.shape[2]),axis=0) for MU, data in neigh_data.iteritems()}
    min_,max_ = 4,9
    fig, ax = plt.subplots()
    index = np.arange(max_+1-min_)
    bar_width=0.7/len(MU_vals)
    rects = [ax.bar(index+bar_width*(i-1),neigh_data[MU][min_:max_+1]/100., bar_width,label=r'$\mu=%.1f$'%-MU) for i,MU in enumerate(MU_vals)]
    plt.legend(loc='best',frameon=False)
    ax.set_xticklabels(np.arange(min_-1,max_+1))
    formatting('Polygon side #','Frequency')
    plt.savefig('side_distr.pdf')

# plot_force(MU_vals)
plot_neigh_distr(MU_vals)
