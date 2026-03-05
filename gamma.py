# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 07:59:56 2016

@author: 331884
"""

import numpy as np
from matplotlib import pyplot as plt
#import scipy as sp
from scipy.interpolate import pchip
def gamma(ref_x, ref_prof,meas_x,meas_prof,dd,dta,norm,thres,interpinterval):

    local_fit_px=1;
    #interpinterval=0.01;
    # comment out line below for global
    #dd = dd*max(ref_prof);
    #meas_x=ref_x
    #np.where(ref_x==0)     
    
    # added thres to account for thresholding.  Will cut profile to what ever is less than thres prior to norm
    #Comment out np.delete lines to do it post norm
    ref_prof = pchip(ref_x, ref_prof)(np.arange(min(ref_x),max(ref_x),interpinterval));
    ref_x = np.arange(min(ref_x),max(ref_x),interpinterval);
    meas_prof = pchip(meas_x, meas_prof)(np.arange(min(meas_x),max(meas_x),interpinterval));
    meas_x = np.arange(min(meas_x),max(meas_x),interpinterval)
    #make profiles same lengyj
    meas_prof= np.delete(meas_prof,np.where(meas_x<min(ref_x)))
    meas_x = np.delete(meas_x,np.where(meas_x<min(ref_x)))
    meas_prof= np.delete(meas_prof,np.where(meas_x>max(ref_x)))
    meas_x = np.delete(meas_x,np.where(meas_x>max(ref_x)))
    ref_prof= np.delete(ref_prof,np.where(ref_x<min(meas_x)))
    ref_x = np.delete(ref_x,np.where(ref_x<min(meas_x)))
    ref_prof= np.delete(ref_prof,np.where(ref_x>max(meas_x)))
    ref_x = np.delete(ref_x,np.where(ref_x>max(meas_x)))
   #threshold pre norm
    #meas_x=np.delete(meas_x,np.where(meas_prof<thres))
    #meas_prof=np.delete(meas_prof,np.where(meas_prof<thres))
    #ref_x=np.delete(ref_x,np.where(ref_prof<thres))
    #ref_prof=np.delete(ref_prof,np.where(ref_prof<thres))
    if norm == 1:#norm to CAX
        i1=(np.abs(ref_x)).argmin();i2=(np.abs(meas_x)).argmin() #indices for normilzation
        ref_prof = ref_prof/(np.mean(ref_prof[i1-2:i1+2]));
        meas_prof= meas_prof/(np.mean(meas_prof[i2-2:i2+2]))#normalize to centeral 5 pixels
    if norm == 2: #norm to dmax
        ref_prof = ref_prof/(np.max(ref_prof)); 
        meas_prof= meas_prof/(np.max(meas_prof));
    if norm==3: #norm to global max
        gmax = max(max(meas_prof),max(ref_prof));
        ref_prof = ref_prof/gmax
        meas_prof= meas_prof/gmax;
    #threshold post norm
    ref_x=np.delete(ref_x,np.where(ref_prof<thres))
    ref_prof=np.delete(ref_prof,np.where(ref_prof<thres))
    meas_x=np.delete(meas_x,np.where(meas_prof<thres))
    meas_prof=np.delete(meas_prof,np.where(meas_prof<thres)) 
    dd_ref = np.reshape(ref_prof/dd,(len(ref_prof),1));
    dta_ref = np.reshape(ref_x/dta,(len(ref_x),1));
    dd_meas = np.reshape(meas_prof/dd,(len(meas_prof),1));
    dta_meas = np.reshape(meas_x/dta,(len(meas_x),1));
    refL=len(dta_ref);
    measL=len(dta_meas)
    
    if local_fit_px <=1:
        #RY = np.transpose(np.multiply(dd_ref,np.transpose(np.ones([1,measL]))));
        #MY = np.ones([refL,1])*np.transpose(dd_meas)
        #MY = (np.ones([refL,1])*np.transpose(dd_meas));
        #RX = np.transpose(np.multiply(np.transpose(np.ones([1,measL])),dta_ref));
        #MX = np.ones([refL,1])*np.transpose(dta_meas);
       # Debugging output to ensure consistent values
        #print(f"refL: {refL}, measL: {measL}")
        #print(f"dta_ref shape: {dta_ref.shape}, dta_meas shape: {dta_meas.shape}")
        
        # Adjust ones_matrix to have correct shapes that match dta_ref and dta_meas
        ones_matrix = np.ones((refL, measL))  # Correct dimensions based on refL and measL
        
        # Perform element-wise multiplications
        RY = np.multiply(dd_ref, ones_matrix)  # (refL, measL)
        MY = np.multiply(np.ones((refL, 1)), dd_meas.T)  # Adjust dd_meas to align with refL
        RX = np.multiply(dta_ref, ones_matrix)  # Ensure correct shape matching with dta_ref
        MX = np.multiply(np.ones((refL, 1)), dta_meas.T)  # Correct alignment of dta_meas


        
    gint=np.sqrt(np.square(RY-MY) + np.square(RX-MX));
    Gammas=np.min(gint,axis=0)
    #print(MY)
    #Gammas = np.min(np.sqrt(np.power((RY-MY),2) +np.power((RX-MX),2)),axis=0);
    #if len(_x)==len(Gammas):
    if len(meas_x) <= len( ref_x):   
       gamma_x = meas_x
    else:
        gamma_x = ref_x  
        
    if len(gamma_x) != len(Gammas):
        idif=int(np.abs(len(gamma_x)-len(Gammas))/2)
        idex1=np.arange(0,idif,1)  
        Gammas=np.delete(Gammas,idex1)
        idex2=np.arange(len(Gammas)-idif,len(Gammas),1)
        Gammas= np.delete(Gammas,idex2)
        #print(str(len(gamma_x)) +  str(len(Gammas))) 
    if len(gamma_x) < len(Gammas):
         Gammas=np.delete(Gammas,0)
    if len(gamma_x) > len(Gammas):
        gamma_x = np.delete(gamma_x,0)
         
            #gamma_x(dta_ref>max(dta_meas))=[];
            #Gammas(dta_ref>max(dta_meas))=[];
            #gamma_x(dta_ref<min(dta_meas))=[];
            #Gammas(dta_ref<min(dta_meas))=[];
    #        return Gammas, gamma_x
    #(Gammas, gamma_x)=profilegamma(ref_y,x,meas_y,x,dta,dd,1)
    #display results
    #mg=np.mean(Gammas);maxg=np.max(Gammas);stdg=np.std(Gammas)
    #gamma_x=meas_x;
   # print('Points: ' + str(len(Gammas)-np.count_nonzero(np.isnan(Gammas))))
    #print(['Passing Rate: ' + str(sum(Gammas<=1)/(len(Gammas)-np.count_nonzero(np.isnan(Gammas))))])
    '''
    plt.figure(3)
    plt.subplot(211)
    plt.title("Profiles")
    plt.plot(ref_x,ref_prof,'.r')
    plt.plot(meas_x,meas_prof,'.k')
    plt.ylabel('Relative Dose')
    plt.xlabel('Off axis distance [mm]')
    plt.ylim(0,np.max(ref_prof)+.1)
    plt.subplot(212)
    plt.plot(gamma_x,Gammas,'.k')
    plt.ylabel('Gamma')
    plt.xlabel('Off Axis Distance [mm]')
    plt.ylim(np.min(Gammas)-0.10,np.max(Gammas)+1)
    plt.show()'''
    return (gamma_x,Gammas)