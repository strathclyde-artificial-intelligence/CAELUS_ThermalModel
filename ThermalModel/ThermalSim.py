import numpy as np
from scipy.integrate import odeint
from .ODEMatrices import ODEMatrices

SIM_STATE_SPACE_SIZE = 6

class ThermalSim():
    
    def __init__(self, params, atmos_model, thermal_model):
        self.__params = params
        self.__atmos_model = atmos_model
        self.__thermal_model = thermal_model
        self.ode_matrices = ODEMatrices(params)

    def solve(self,
        time_start,
        time_end,
        max_step = 1,
        initial_state=[0]*SIM_STATE_SPACE_SIZE):

        t = np.linspace(time_start, time_end)
        sol = odeint(self.__thermal_model, initial_state, t, args=(self.ode_matrices, self.__atmos_model, self.__params), tfirst=True, hmax=max_step)
        return t, sol