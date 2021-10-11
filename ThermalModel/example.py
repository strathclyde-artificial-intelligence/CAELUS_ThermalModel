# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
from wrapper import container_wrapper



# RUN EXAMPLE

T0 = [0,0,0,0,0,0]
time_end = 8000.0

t, sol = container_wrapper(T0, time_end)
    
    
# plot results

plt.plot(t,sol)
plt.xlabel('time')
plt.ylabel('y(t)')
plt.show()