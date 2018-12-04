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
    
def plot_graph_type(filename,ax,color,marker='o',label=None,error=False,bestfit=False,start_index=None,end_index=None,bvals_to_use=None,alpha=1.):
    fix_data = np.loadtxt(filename,skiprows=1,dtype=float,comments='#')
    fix_data[:,1]= fix_data[:,1]/(fix_data[:,2]+fix_data[:,1])
    fix_data = fix_data[:,:2].T
    b_vals = fix_data[0]
    fixprobs=fix_data[1]
    if bvals_to_use is not None: 
        mask = np.array([np.where(b_vals==b)[0][0] for b in bvals_to_use])
        b_vals,fixprobs = b_vals[mask],fixprobs[mask]
    ax.plot(b_vals[start_index:end_index],fixprobs[start_index:end_index]*100,linestyle='',marker=marker,color=color,label=label,alpha=alpha)
    if error:
        confidence_intervals_static = [confint(p,10000)*100 for p in fix_data[1]]
        ax.errorbar(fix_data[0][start_index:end_index],fix_data[1][start_index:end_index]*100,yerr=confidence_intervals_static[:7],fmt='none',elinewidth=1.5,color=color)
    if bestfit:
        slope,intercept,r,p,st = linregress(fix_data[0][start_index:end_index],fix_data[1][start_index:end_index])
        ax.plot(fix_data[0][start_index:end_index],slope*fix_data[0][start_index:end_index]*100+intercept*100,ls='--',color=color)
    
def plot_vt_type(folder,ax,color,marker='o',label=None,error=False,bestfit=False,bvals_to_use=None,bestfit_bvals=None):
    filenames = [filename for filename in os.listdir(folder)if filename[:3]=='fix']
    results = np.array([read_data('%s/%s'%(folder,filename)) for filename in filenames])
    b_vals = np.array([filename[3:] for filename in filenames],dtype=float)
    fixprobs = results[:,0]/results[:,1]
    if bvals_to_use is not None: 
        mask = np.array([np.where(b_vals==b)[0][0] for b in bvals_to_use])
        b_vals,fixprobs = b_vals[mask],fixprobs[mask]
    ax.plot(b_vals,fixprobs*100,marker=marker,ls='',linewidth=1.5,label=label,color=color)
    if error:
        confidence_intervals = [confint(p,r.sum())*100 for p,r in zip(fixprobs,results)]
        ax.errorbar(b_vals,fixprobs*100,yerr=confidence_intervals,fmt='none',elinewidth=1.5,color=color)
    if bestfit:
        if bestfit_bvals is not None:
            mask = np.array([np.where(b_vals==b)[0][0] for b in bestfit_bvals])
            slope,intercept,r,p,st = linregress(b_vals[mask],fixprobs[mask])
        else: slope,intercept,r,p,st = linregress(b_vals,fixprobs)
        ax.plot(b_vals,slope*b_vals*100+intercept*100,ls='--',color=color)
        print 'critical ratio for '+ folder+':'
        print (0.01-intercept)/slope

def plot_theory(filename,ax,color,start_index=None,end_index=None,label=None,percent=False,zorder=None):
    fix_theory = np.loadtxt(filename).T
    if not percent: ax.plot(fix_theory[0][start_index:end_index],fix_theory[1][start_index:end_index]*100,color=color,label=label,zorder=zorder)
    else: ax.plot(fix_theory[0][start_index:end_index],fix_theory[1][start_index:end_index],color=color,label=label,zorder=zorder)

sns.set_style('white')    

cp = sns.color_palette('deep6')

#Figure 1 
fig,ax = plt.subplots()
plot_neutral_line(ax,2,12)
plot_theory('EGTpd_av_db/allen_result/vt_graph_6',ax,cp[0],end_index=None)
plot_theory('EGTpd_av_db/allen_result/hex_graph',ax,cp[1],end_index=None)
plot_graph_type('EGTpd_av_db/no_migration/vt_graph_6',ax,cp[0],marker='s',error=False,end_index=None,alpha=1.0,label='DT')
plot_graph_type('EGTpd_av_db/no_migration/hex_graph',ax,cp[1],marker='o',error=False,end_index=None,alpha=1.0,label='HL')
plt.yticks(np.arange(0.0,2.5,0.5))
formatting(r'Benefit-to-cost ratio, $b/c$',r'Fixation Probability, $\rho _C$ (%)',large=False,legend=True)
plt.savefig('hex_vg6.pdf')

#Figure 5
fig,ax = plt.subplots()
plot_neutral_line(ax,1.5,4.0)
plot_theory('VTpd_av_decoupled_theory',ax,cp[0],percent=True,end_index=6)
plot_vt_type('VTpd_av_decoupled',ax,cp[0],error=False,marker='o',bestfit=False,bvals_to_use=np.arange(1.5,4.25,0.25))
plt.yticks(np.arange(0.4,2.0,0.2))
formatting(r'Benefit-to-cost ratio, $b/c$',r'Fixation Probability, $\rho _C$ (%)',large=False,legend=True)
plt.savefig('vtcompare.pdf')

#Figure 6
fig,ax = plt.subplots()
plot_neutral_line(ax,2,8)
plot_theory('EGTpd_av_db/allen_result/hex_graph',ax,cp[1],end_index=7,label='HL (death-birth)')
plot_vt_type('VTpd_av_decoupled',ax,cp[0],error=False,marker='o',bestfit=False,label='VT (decoupled)',bvals_to_use=[2,2.5,3,3.5,4,5,6,8])
plot_theory('VTpd_av_decoupled_theory',ax,cp[0],percent=True,start_index=1)
plot_vt_type('VTpd_av_db',ax,cp[2],error=False,marker='s',bestfit=True,label='VT (death-birth)',bestfit_bvals=[6.,6.5,7,7.5,8])
plt.yticks(np.arange(0.0,4.5,1.0))
formatting(r'Benefit-to-cost ratio, $b/c$',r'Fixation Probability, $\rho _C$ (%)',large=False,legend=True)
plt.savefig('vtcompare.pdf')

plt.show()

