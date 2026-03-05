# -*- coding: utf-8 -*-
"""
Created on Sat Nov  3 09:29:26 2018
calc dose difference of interpolated profiles 
For dose difference it normalizes to dmax of each and subtracts
@author: nknutson
"""
from scipy import interpolate as interp
from matplotlib import pyplot as plt
import numpy as np
def dosedif(x1, d1, x2, d2, norm=2):
    if norm == 1:  # Normalize to CAX
        i1 = (np.abs(x1)).argmin()
        i2 = (np.abs(x2)).argmin()
        d1 = d1 / (np.mean(d1[i1-2:i1+2]))
        d2 = d2 / (np.mean(d2[i2-2:i2+2]))
    elif norm == 2:  # Normalize to Dmax of each profile
        d1 = d1 / np.max(d1)
        d2 = d2 / np.max(d2)
    elif norm == 3:    #norm to global max
        gmax = np.max([np.max(d1), np.max(d2)])

        d1=d1/gmax
        d2=d2/gmax
    d1int = interp.pchip(x1, d1)    
    d2int = interp.pchip(x2, d2)
    xint = np.arange(max(min(x1), min(x2)), min(max(x1), max(x2)), 0.01)
    dd = d1int(xint) - d2int(xint)
   
        
    return (xint, dd)
def dta2(x1,d1,x2,d2):#old slow method Dont use any more use Dta below that is 10 x faster
    ddl=[];sl=[];dtal=[];dtaval=[]
    d1int=interp.pchip(x1,d1)    
    d2int=interp.pchip(x2,d2)
    stest=np.arange(-.35,.35,0.01)
    xint=np.arange(max(min(x1),min(x2)),min(max(x1),max(x2)),0.01)
    for x in xint:
        dtal=[];sl=[];ddl=[];
        for s in stest:
            #print(s)
            ddpoint=np.abs(d1int(x)-d2int(x+s))
            ddl.append(ddpoint)
            sl.append(s)
        index=np.argmin(ddl)
        dtal.append(np.abs(sl[index]))
        dtaval.append(min(dtal))
        #print(dtaval)
        #dtal=np.abs(dtal)
    dtaval=np.asarray(dtaval)
    return (xint,dtaval)   

def dta(x1, d1, x2, d2,dta):
    # same interpolators you already use
    d1int = interp.pchip(x1, d1)
    d2int = interp.pchip(x2, d2)

    # same grids you already use
    xint  = np.arange(max(min(x1), min(x2)),
                      min(max(x1), max(x2)), 0.01)
    stest = np.arange(-(dta+0.1), dta+0.11, 0.01)  #test the Dta threshold plus 1mm

    # evaluate once on vectorized grids (no Python loops)
    d1_on_x    = d1int(xint)                                  # (N,)
    d2_on_grid = d2int(xint[:, None] + stest[None, :])        # (N, S)

    # find shift with minimum |dose diff| for each x
    diffs   = np.abs(d1_on_x[:, None] - d2_on_grid)           # (N, S)
    best_ix = np.argmin(diffs, axis=1)                        # (N,)
    dtaval  = np.abs(stest[best_ix])                          # (N,) in cm

    return xint, dtaval

'''#test profiles
x1=np.arange(-2.1,2.3,.25)
x2=np.arange(-2,2.4,.01)
d1=np.exp(-np.power(x1 , 2.) / (2 * np.power(.5, 2.)))
d2=np.exp(-np.power(x2 - .2, 2.) / (2 * np.power(.5, 2.)))
n1 = np.random.normal(0,.005,len(d1))
n2 = np.random.normal(0,.005,len(d2))
d1=d1+n1;d2=d2+n2;
testx,testy=dosedif(x1,d1,x2,d2)
plt.plot(x1,d1,'.b');plt.plot(x2,d2,'.g');plt.plot(testx,testy,'.k')
'''