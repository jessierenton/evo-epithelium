import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sns.set_style('white')

def save_averaged_coop_interaction(data_CC):
    coop_int = np.array([np.mean(dat_n)/(i+1.) for i,dat_n in enumerate(data_CC)])
    np.savetxt('sigma_cc',coop_int)

def load_data(folder,dtype=float):    
    return [np.loadtxt(folder+'/n_%d'%n,dtype=dtype) for n in range(1,N)]
    
def increase_prob(b,n,data_CC_n,data_CD_n):
    return (1.-float(n)/N)*float(n)/N*(1+DELTA*av_coop_payoff(b,n,data_CC_n))/(1+DELTA*av_payoff(b,n,data_CC_n,data_CD_n))

def decrease_prob(b,n,data_CC_n,data_CD_n):
    return float(n)/float(N)*(1-float(n)/N*(1+DELTA*av_coop_payoff(b,n,data_CC_n))/(1+DELTA*av_payoff(b,n,data_CC_n,data_CD_n)))

def av_coop_payoff(b,n,data_CC_n):
    return -c +b*data_CC_n/n

def av_payoff(b,n,data_CC_n,data_CD_n):
    return (-c*n+b*(data_CC_n+data_CD_n))/float(N)

def gos(b,data_CC,data_CD):
    return [increase_prob(b,i+1,data_CC_n,data_CD_n)-(b,i+1,data_CC_n,data_CD_n) for i,(data_CC_n,data_CD_n) in enumerate(zip(data_CC,data_CD))]

def gos_weak_selection(b,sigma_cc):
    n=np.arange(1,100,dtype=float)
    return np.hstack([0,(n/N**2)*DELTA*(-c*(N-n)+b*(sigma_cc*N-n)),0])

def plot_gradient_mean_std(gradient_data,ax,show_error=True,label=None):
    x = np.arange(1,N)
    gradient_mean = np.array([np.mean(grad_n) for grad_n in gradient_data])
    ax.plot(x/100.,gradient_mean,label=label)
    if show_error:
        gradient_std = np.array([np.std(grad_n) for grad_n in gradient_data])
        ax.fill_between(x/100.,gradient_mean-gradient_std,gradient_mean+gradient_std,alpha=0.3)

def save_gradient_means(gradient_data):
    np.savetxt('av_gos',[np.mean(grad_n) for grad_n in gradient_data])
    
def cube_fit(gradient_data,ax,label=None,fit_region=(None,None)):
    x = np.arange(1,N)/100.
    gradient_mean = np.array([np.mean(grad_n) for grad_n in gradient_data])
    gradient_std = np.array([np.std(grad_n) for grad_n in gradient_data])
    ax.plot(x,gradient_mean,label=label)
    fit_params = np.polyfit(x[fit_region[0]:fit_region[1]],gradient_mean[fit_region[0]:fit_region[1]],deg=3)
    poly = np.poly1d(fit_params)
    ax.plot(x,poly(x))
    return poly

def theory_plot(S,T):
    x = np.arange(1,N)/100.
    plt.plot(x,DELTA*x*(1-x)*(x*(1-S-T)+S))
    
def fix_prob(b,c,sigma_cc):
    f = lambda j: -c/N +b/N*(sigma_cc[j-1]*N-j)/(N-j)
    return 1./N +DELTA/N*sum(sum(f(j) for j in range(1,k+1)) for k in range(1,N))

def formatting(xlabel,ylabel,large=False):
    if large: labelsize,ticksize = 26,18
    else: labelsize,ticksize = 16,12   
    plt.xlabel(xlabel,size=labelsize,labelpad=10)
    plt.ylabel(ylabel,size=labelsize,labelpad=10)
    plt.xticks(size=ticksize)
    plt.yticks(size=ticksize)
    plt.tight_layout()

N = 100
DELTA = 0.025
c = 1.
if __name__ == '__main__':    
    
    sigma_cc = np.loadtxt('sigma_cc')
    x = np.arange(0.01,1.0,0.01)
    sns.set_context('paper')
    fig,ax = plt.subplots()
    for b in np.arange(1.5,4.5,0.5):
        plt.scatter(np.hstack([0,x,1]),gos_weak_selection(b,sigma_cc)*10**(3),label=r'$b=%.1f$'%b,marker='o',s=8)

    plt.plot(np.arange(0,1,0.01),[0]*100,ls='--',color='grey')
    plt.legend(loc='best')
    formatting(r'Fraction of cooperators, $n/N$',r'$G^A(n)$ $(x10^{-3})$')
    plt.savefig('gos.pdf')
    plt.show()
    plt.plot(x,sigma_cc)
    
    
    #
# gradient_data = [increase_prob(b,i+1,data_CC[i],data_CD[i])-decrease_prob(i+1) for i in range(0,N-1)]
# save_gradient_means(gradient_data)

# for b in (1.5,):
#     gradient_data = [increase_prob(b,i+1,data_CC[i],data_CD[i])-decrease_prob(b,i+1,data_CC[i],data_CD[i]) for i in range(0,N-1)]
#     plot_gradient_mean_std(gradient_data,ax,show_error=False,label=r'$b=%.1f$'%b)
#     theory_plot(-0.0081607199999999998,1.0097669499999999)
#     # poly_fit = cube_fit(gradient_data,ax)
