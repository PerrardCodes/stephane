

import numpy as np


def pressure(h,Dtheta,R,R0=1):
    """
    Return the pressure difference for a spherical cap
    INPUT 
    -----
    h : thickness in inches
    R : Radius
    R0 : initial disk radius
    """
    
    h = h*0.0254 #thickness in meter from h in inches
  #  Dtheta = Dtheta*np.pi/180
    
    E = 5.5*10**6 #Young modulus
    nu = 0.5 #
    mu = E/2/(1+nu)
    l = Dtheta*R/(2*R0)
    
    dP = 2*h/R*mu*(l**2-l**(-4))*1000 #conversion of R from mm to m
    
    return dP/10**5 #return the result in bar
    
    
def main():
    
#experiments from the 12/19
    fx=0.101


    h =0.06
    Dtheta=(360-131)*np.pi/180
    R = 417*fx
    R0 = 376*fx
    dP1 = pressure(h,Dtheta,R,R0=R0)

    h =0.1
    R = 433*fx
    R0 = 376*fx
    Dtheta=2*np.pi-2*np.arcsin(R0/R)
    dP2 = pressure(h,Dtheta,R,R0=R0)

    h =0.2
    R = 815.*fx/2
    R0 = 748.*fx/2
    Dtheta=2*np.pi-2*np.arcsin(R0/R)

    dP3 = pressure(h,Dtheta,R,R0=R0)    
#    print(np.arcsin(R0/R))
    print(dP1,dP2,dP3)
    
#experiment from the 12/26
    fx=0.1863
    
    h = 0.08
    R = 668.*fx/2
    R0 = 409.*fx/2
    
    Dtheta=2*np.pi-2*np.arcsin(R0/R)
    dP4 = pressure(h,Dtheta,R,R0=R0)
    
    print(dP4)
    
#main()