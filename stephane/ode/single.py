import numpy as np
import os
import transport
import pickle

def compute():
    R={}
    Vlist = np.linspace(0,6,6*2+1)

    for vr0 in Vlist:
        v0 = vr0/U0
    
        thetalist = theta*np.ones(n)
        vlist = v0*np.ones(n)

        for i,ell in enumerate(llist):
            print(theta,v0,ell)
            R[(theta,v0,ell)]={}
            Fr = U0/np.sqrt(g*ell)
    
            tref,xref,zref = transport.traj_dicho(theta,v0,Fr,Fr0,0,Dt,N=N)
            Xlist,Zlist,X,Z = transport.traj_beam_dicho(thetalist,vlist,Fr,Fr0,sigma,Dt,N=N)
        
            for key in ['v0','theta','Dt','L','sigma','Fr0','Fr','N','ell','Xlist','Zlist','n']:#X,Z
                if key in globals().keys():
                    R[(theta,v0,ell)][key]=globals()[key]
                else:
                    R[(theta,v0,ell)][key]=locals()[key]
    return R

def run(savefolder):
    global theta,N,Dt,n,U0,llist,g,Fr0,sigma,L

    g = 9.81
    H = 0.05
    L = 1.6

    N = 10**3
    n = 100

    Ulist = [2.38,4.76,6.67]
    U0 = Ulist[0]
    Fr0 = U0/np.sqrt(g*H)

    theta = np.pi/8
    sigma = 0.1

    T = L/U0
    Dt = 0.1*T

    nl=3
    llist = np.arange(0.02,1,nl)
    p={}

    R = compute()

    filename = savefolder+'test.pkl'
    with open(filename, 'wb') as f1:
        pickle.dump(R, f1)


def main():
    import time
    savefolder = '/Users/stephane/Documents/Programming/Python/Notebooks/Drop_Turbulence/Mar23/Data/'
    if not os.path.exists(savefolder):
        os.makedirs(savefolder)
    t1 = time.time()
    run(savefolder)
    t2 = time.time()
    print(str(np.round(t2-t1,decimals=1))+' s')
    print('done')

main()
    
