#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  9 16:40:24 2021

@author: gianlucafilippi
"""



import numpy as np

pi = 3.14











"""
    SpreadingResistance(state)

"""
def SpreadingResistance(r1, r2, t, k, h):

    ϵ = r1/r2
    τ = t / r2
    Bi = (h * r2) / k
    λ = pi + 1/(1/ (ϵ * pow(pi,0.5)) )
    ϕ = (np.tanh(λ * τ) + λ/Bi) / (1+ λ/Bi * np.tanh(λ*τ))
    ψ = (ϵ*τ)/pow(pi,0.5) + 1/pow(pi,0.5)*(1-ϵ)*ϕ

    Rsp = ψ / (k*r1*pow(pi,0.5))
    return Rsp










"""
    CalculateMatrices(state)
sfab fdbabe
"""
def CalculateMatrices(dd):


    # ------------- initialisation ------------------------
    # matrix initialisation
    # H = zeros(6,6)  # convection
    # K = zeros(6,6)  # conduction
    # R = zeros(6,6)  # radiation
    # Q = zeros(6)    # external heat source

    # # M = zeros(5)    # vector of masses
    # C = zeros(6)    # vector of heat capacity
    # -----------------------------------------------------





    # -------------- heat capcity ------------------------
    # node 1 : external wall / payload
    C1 = dd["container_V"] * dd["container_rho"] * dd["container_cp"]   # J K-1

    # node 2 : internal wall / payload
    C2 = dd["container_V"] * dd["container_rho"] * dd["container_cp"]   # J K-1

    # node 3 : payload
    C3 = dd["payload_V"] * dd["payload_rho"] * dd["payload_cp"]   # J K-1

    # node 4 : PCM
    C4 = 0
    
    # node 5 : air
    C5 = dd["air_V"] * dd["air_rho"] * dd["air_cp"]   # J K-1
    
    C = np.array([C1, C2, C3, C4, C5, 0])
    # -------------------------------------------------

    Dx_container = dd["Dx_container"]



    # -------- spreading resistances--------

    # from "square" to "circular" geometry
    r_contianer = pow(dd["container_Aout"] / pi, 0.5)
    r_brick = pow(dd["payload_A_cond"] / pi, 0.5)
    r_PCM = pow(dd["PCM_A_cond"] / pi, 0.5)

    R_c2b = SpreadingResistance(r_brick, r_contianer, Dx_container, dd["container_k"], dd["air_hf"])
    R_PCM2b = SpreadingResistance(r_brick, r_PCM, dd["PCM_z"], dd["PCM_ks"], dd["air_hn"])

    # R_c2b=0
    # R_PCM2b=0
    # ------------------------------------------


    # ----------------- conductivity ------------------------

    # node 1 (external wall) <-> node 2 (internal wall/payload)
    K_1_2 = dd["container_k"] * dd["container_Aout"] / Dx_container
    K_2_1 = K_1_2

    # node 2 (internal wall/payload)  <-> node 3 (payload)
    R_23 = pow( dd["payload_k"] * dd["payload_A_cond"] / (dd["payload_z"]/2) , -1)
    R_tot = R_23+R_c2b
    K_tot = pow(R_tot, -1)
    
    K_3_2 = K_tot
    K_2_3 = K_3_2 # R_23^-1

    # node 3 (payload)  <-> node 4 (PCM)                 N.B.! ONLY PAYLOAD PART
    R_34 = pow( dd["payload_k"] * dd["PCM_A_cond"] / (dd["payload_z"]/2),  -1)
    R_tot = R_34+R_PCM2b
    
    K_4_3 = pow(R_34, -1)
    K_3_4 = K_4_3 # R_tot^(-1)
    
    Kpartial = K_4_3
    
    
    K = np.array([ [0, K_1_2, 0, 0, 0, 0], [K_2_1, 0, K_2_3, 0, 0, 0], [0, K_3_2, 0, K_3_4, 0, 0], [0, 0, K_4_3, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0] ])
    # ---------------------------------------------------------



    # --------- convection -------------
    # node 1 (internal wall)  <-> node 6 (external air)
    H_1_6 = dd["air_hn"] * dd["container_Aout"]
    H_6_1 = H_1_6

    # node 2 (internal wall)  <-> node 5 (air)
    H_2_5 = dd["air_hn"] * dd["container_Ain_conv"]
    H_5_2 = H_2_5

    # node 3 (payload)  <-> node 5 (air)
    H_3_5 = dd["air_hn"] * dd["payload_A_conv"]
    H_5_3 = H_3_5

    # node 4 (PCM)  <-> node 5 (air)
    H_4_5 = dd["air_hn"] * dd["PCM_A_conv"]
    H_5_4 = H_4_5
    
    Hpartial = H_4_5
    
    H = np.array([[0, 0, 0, 0, 0, H_1_6], [0, 0, 0, 0, H_2_5, 0], [0, 0, 0, 0, H_3_5, 0], [0, 0, 0, 0, H_4_5, 0], [0, H_5_2, H_5_3, H_5_4, 0, 0], [H_6_1, 0, 0, 0, 0, 0] ])
    # ----------------------------------



    # radiation
    # ----------------------------------
    R = np.array([ [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0] ])
    #  ----------------------------------
     
    
    Q = np.array([0,0,0,0,0])
    
    
    
    
    
    out = {
      "C": C,  
      "K": K,  
      "H": H,  
      "R": R,  
      "Q": Q,  
      "Kpartial": Kpartial,
      "Hpartial": Hpartial
      }
    
    
    return out


# ------------------------------------------------------------------------------