import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 
from shapely.geometry import LineString


#L = [8, 16, 32, 64]
L = [16, 24, 32, 40, 64]
#L=[8]

temperatures = [0.4556158469036384, 0.453, 0.4157943265438966, 0.422, 0.4927759800470331] 
rainbow_colors = ['red', 'orange', 'deepskyblue', 'green', 'blue', 'indigo', 'violet', 'brown']

i=0 

for l in L: 
    
    plt.plot ( l , temperatures[i], 'ro', label=f'L={l}', color=rainbow_colors[i % len(rainbow_colors)])
    plt.plot ()

    i += 1



plt.grid(True)
plt.legend()
plt.title('Temperature of Z2 transition' )
plt.xlabel(r'$\frac{1}{L}$')
plt.ylabel(r'$T_C^{Z2}$')

plt.show()