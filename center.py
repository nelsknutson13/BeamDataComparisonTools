# -*- coding: utf-8 -*-
"""
Created on Thu Nov  1 09:33:16 2018
PAss profile and the threshold to center around 0.5 = 50%
@author: nknutson
"""
import numpy as np
from scipy import interpolate as interp
import warnings

def center(x1,d1,th):
     dtest=interp.pchip(np.asarray(x1),np.asarray(d1))
     xinterp=np.arange(np.min(x1),np.max(x1),((x1.iloc[1]-x1.iloc[0])/10))
     dtest=dtest(xinterp)
     max_d = max(dtest)  # Find the maximum y value
     xs = [x for x in range(len(dtest)) if dtest[x] > max_d*th]# Find pos of half max on each side
     if not xs:
        warnings.warn(f"No points found above threshold {th}. Returning 0.")
        return 0
     x1=xinterp[min(xs)]#LEft Side max
     x2=xinterp[max(xs)]#Right Side Max
     s1=((np.abs(x1)-np.abs(x2)))/2# shift to make them equal
     return(s1)
    