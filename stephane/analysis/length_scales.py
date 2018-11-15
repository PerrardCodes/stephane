
import numpy as np

def dissipation_rate(Omega,d=1):
    factor = 3./d
    
    nu = 1  # in mm^2/s
    Enstrophy = factor*Omega**2
    result = nu * Enstrophy
    print('Dissipation rate : '+str(avg(result,3)) + '  mm^2/s^-3')
    return result
    
def K_scale(epsilon):
    nu = 1  # in mm^2/s
    result = (nu**3/epsilon)**0.25
    
    print('Kolmogorov scale : '+str(avg(result,3)) + '  mm')
    return result
    
def T_scale(Omega,E):
    E_rms = np.nanmean(E)
    Enstrophy = np.nanmean(Omega**2)
    
    result = np.sqrt(1./2.*E_rms/Enstrophy)
    
    print('Taylor scale : '+str(avg(result,3)) + '  mm')
    return result
    
def Re(micro,eta):
    return (micro/eta)**4  # wrong in general !!
    
def Re_lambda(E,micro):
    nu = 1 # in mm^2/s
    return np.sqrt(E)*micro/nu
    
def I_scale(Re,E):
    nu = 1 # in mm^2/s
    return nu*Re/E
    
def vs_t(M,fun,Dt=1):
    ny,nx,nt = M.shape()
    
    for i in range(Dt,nt-Dt):
        S[i] = fun()        
    
def avg(a,n):
    return int(a*10**n)*1./10**n
    
def main():
    Dt = 10
    tax = range(1,1500,1)

    E_mean = np.zeros(len(tax))
    dU_mean = np.zeros(len(tax))
    dU_std = np.zeros(len(tax))

    figs = {}
    display_part = False
    for i,t in enumerate(tax):
        strain = M_2015_09_17[0].strain[...,t:t+Dt]
        E = M_2015_09_17[0].E[...,t:t+Dt]
    
        E_mean[i] = np.nanmean(E)
        dU_mean[i] = np.nanmean(strain)
        dU_std[i] = np.nanstd(strain)
    
        if display_part:
            xbin,n = graphes.hist(strain,fignum=1)    
            figs.update(graphes.legende('dU/dx','PDF',''))
    
            graphes.graph(xbin/dU_mean[i],n/max(n),fignum=2)
            figs.update(graphes.legende('dU/dx','Normalized PDF',''))
            graphes.save_graphes(M_2015_09_17[0],figs,prefix='Strain/',suffix='Dt_'+str(Dt)+'_'+str(t))
        #.shape
    tf = np.asarray(M_2015_09_17[0].t)[tax]-0.04

    U_rms = np.sqrt(E_mean)
    l_micro = np.sqrt(E_mean/(3/2*dU_mean)**2)
    Re_l = U_rms * l_micro

    graphes.graphloglog(tf,l_micro,label='ko',fignum=3)
    figs.update(graphes.legende('t (s)','lambda',''))

    graphes.graphloglog(tf,U_rms,label='ro',fignum=4)
    figs.update(graphes.legende('t (s)','U_rms',''))


    graphes.graphloglog(tf,Re_l,label='k+',fignum=5)
    figs.update(graphes.legende('t (s)','Re_lambda',''))

    #graphes.graphloglog(tf,100*tf**(-3/2),label='r',fignum=3)
    #graphes.set_axis(10**-2,10**2,10**-2,10**4)


    print(dU_std/dU_mean)
    graphes.semilogx(np.asarray(M_2015_09_17[0].t)[tax]-0.04,dU_std/dU_mean,label='k^',fignum=6)
    figs.update(graphes.legende('t (s)','omega>_rms/<omega>',''))
    #graphes.set_axis(10**-2,10**2,0,0.8)

    graphes.save_graphes(M_2015_09_17[0],figs,prefix='Vorticity',suffix='log')