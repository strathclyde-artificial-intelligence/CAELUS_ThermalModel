#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 18:56:09 2021

@author: gianlucafilippi
"""



import numpy as np
from scipy.integrate import odeint





from container_inputs import input_geometry 
from container_helpfunctions import CalculateMatrices 
from container_models import model_atmosferic, container_thermalnetwork 




"""
wrapper function
"""
def container_wrapper(T0, time_end):




    ## define parameters
    inputs = input_geometry()


    # matrices
    M = CalculateMatrices(inputs)



    # time span
    time_0 = 0.0
    # time_end = 8000.0
    # tspan = (time_0, time_end)
    t = np.linspace(time_0, time_end)




    ## solve ODE
    sol = odeint(container_thermalnetwork, T0, t, args=(inputs, M, model_atmosferic, ))



    # define problem
    #prob = ODEProblem(Container_2.container_thermalnetwork!, T0, tspan, params)
    ## solve ODE
    #sol = solve(prob);

    # return sol.u[end]
    return t, sol





# RUN EXAMPLE

# T0 = [0,0,0,0,0,0]
# time_end = 8000.0

# t, sol = container_wrapper(T0, time_end)
    
    
# # plot results

# plt.plot(t,sol)
# plt.xlabel('time')
# plt.ylabel('y(t)')
# plt.show()