#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  9 16:40:24 2021

@author: gianlucafilippi
"""

import numpy as np




"""

atmosferictemperature

"""
def model_atmosferic(t):

    if (t< 12*60):        
        T = 14
    if (12*60<=t< 45*60):
        T = 12
    if (45*60<=t< 105*60):
        T = 16
    if (105*60<=t< 131*60):
        T = 12
    if (t> 131*60):
        T = 16

    return T






"""
MAIN ODE FUNCTION
"""

def container_thermalnetwork(T, t, inputs, M, model_atmosferic):

    
    dT = [0,0,0,0,0,0]
    
    #inputs, M, model_atmosferic = params
    
    C = M["C"]
    K = M["K"]
    H = M["H"]
    R = M["R"]
    Q = M["Q"]
    
    Kpartial = M["Kpartial"]
    Hpartial = M["Hpartial"]

    
    # C, K, H, R, Q, Kpartial, Hpartial, pcm_prop, atmmodel = params


    # ambient temperature from Atmosferic model
    T[5] = model_atmosferic(t)


    # PCM interpolation
    K4 = []
    L = (T[3]-inputs["PCM_t1"]) / (inputs["PCM_t2"]-inputs["PCM_t1"])
    # PCM properties
    if L < 0:
        Ln = 0
        Mass_pcm = inputs["PCM_rho"] * inputs["PCM_V"]
        cp = inputs["PCM_cps"]
        K4 = inputs["PCM_ks"] * inputs["PCM_A_cond"] / (inputs["PCM_z"]/2)

    elif (0 <= L < 1):
        Ln = L
        Mass_pcm = inputs["PCM_rho"] * inputs["PCM_V"]
        cp = (inputs["PCM_cpl"] * (T[3] - inputs["PCM_t1"]) + inputs["PCM_h"] + inputs["PCM_cps"] * ((inputs["PCM_t2"] - T[3])) ) / (inputs["PCM_t2"] - inputs["PCM_t1"])
        # keff =   (pcm_prop.k_l*(T[4]-pcm_prop.Tm1) + pcm_prop.k_s*(pcm_prop.Tm2-T[4]) )/ (pcm_prop.Tm2-pcm_prop.Tm1)
        keff = inputs["PCM_ks"] + L * (inputs["PCM_kl"] - inputs["PCM_ks"])
        K4 = keff * inputs["PCM_A_cond"] / (inputs["PCM_z"]/2)

    else:
        Ln = 1
        Mass_pcm = inputs["PCM_rho"] * inputs["PCM_V"]
        cp = inputs["PCM_cpl"]
        K4 = inputs["PCM_kl"] * inputs["PCM_A_cond"] / (inputs["PCM_z"]/2)

    

    K4new = pow((1/Kpartial + 1/K4), -1)
    K[2,3] = K4new
    K[3,2] = K4new

    C[3] = Mass_pcm * cp

    H45new = pow((1/Hpartial + 1/K4), -1)
    H[3,4] = H45new
    H[4,3] = H45new




    TT_row = np.array([ [T[0]]*6, [T[1]]*6, [T[2]]*6, [T[3]]*6, [T[4]]*6, [T[5]]*6 ])
    TT_col = TT_row.transpose()
    
    TT = TT_row - TT_col
    TTR =  pow(TT_row,4) - pow(TT_col,4)
    
    dT = 1./C * ( sum(np.dot(H, TT)) + sum(np.dot(K, TT)) + sum(np.dot(R, TTR)) + Q ) 
    dT[5] = 0

    return dT








