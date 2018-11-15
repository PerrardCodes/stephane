

import scipy.optimize
import numpy as np

def fit(fun,X,Y):
    try:
        popt,pcurv = scipy.optimize.curve_fit(fun,X,Y)
    except ValueError:
        print("NaN values encountered, fit skipped")
           # input()
        pcurv=[]
        popt=[np.nan]
    except RuntimeError:
        print("Fitting did not converge, arbitrarly chosen to 1")
        pcurv = []
        popt = [np.nan]
        
    return popt,pcurv
    
def parabola(t,a):
    return a*t**2
    
def exp(t,tau):
    return np.exp(t*tau)
    
def gaussian(x,a):
    return np.exp(-a*x**2)