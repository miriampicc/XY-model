import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 
from shapely.geometry import LineString

 
"""temperatures = [ 0.009, 0.01, 0.02, 0.03, 0.05, 0.07, 0.08, 0.09, 0.1, 0.11, 0.13, 0.14, 0.15, 
                0.17, 0.2, 0.22, 0.25, 0.27, 0.3, 0.32, 0.35, 0.37, 0.4, 0.45, 0.5, 0.52, 0.55, 
                0.57, 0.6, 0.62, 0.65, 0.67, 0.7, 0.72, 0.75, 0.77, 0.8, 0.82, 0.85, 0.87, 0.9, 
                0.92, 0.94, 0.95, 0.96, 0.98, 1.0, 1.02, 1.04, 1.05, 1.06, 1.07, 1.09, 1.1, 1.13, 1.15, 1.17, 1.2, 1.22, 1.25, 
                1.27, 1.3, 1.4, 1.45, 1.5, 1.7, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5]   """

"""temperatures =[0.3, 0.32, 0.35, 0.36, 0.37, 0.39, 0.4, 0.405, 0.41, 0.415, 0.42, 0.422, 0.424, 0.426, 0.428, 
               0.43, 0.432, 0.434, 0.436, 0.438, 0.44, 0.442, 0.444, 0.446, 0.448, 0.45, 0.455, 0.46, 0.465, 
               0.47, 0.475, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6,
                0.62, 0.65, 0.67, 0.7, 0.72, 0.75, 0.77, 0.8, 0.82, 0.85, 0.87, 0.9, 0.92, 0.94, 0.95, 0.96, 0.98, 1.0]
"""
"""temperatures = [ 0.009, 0.01, 0.02, 0.03, 0.05, 0.07, 0.08, 0.09, 0.1, 0.11, 0.13, 0.14, 0.15, 
                0.17, 0.2, 0.22, 0.25, 0.27, 0.3, 0.32, 0.35, 0.37, 0.4, 0.45, 0.5, 0.52, 0.55, 
                0.57, 0.6, 0.62, 0.65, 0.67, 0.7, 0.72, 0.75, 0.77, 0.8, 0.82, 0.85, 0.87, 0.9, 
                0.92, 0.94, 0.95, 0.96, 0.98, 1.0, 1.02, 1.04, 1.05, 1.06, 1.07, 1.09, 1.1, 1.13, 1.15, 1.17, 1.2, 1.22, 1.25, 
                1.27, 1.3, 1.4, 1.45, 1.5] """


"""temperatures = [0.4, 0.405, 0.41, 0.415, 0.42, 0.422, 0.424, 0.426, 0.428, 0.43, 0.4305, 0.4310, 0.4315, 0.4320, 0.4325,
                0.4330, 0.4335, 0.4340, 0.4345, 0.4350, 0.4355, 0.4360, 0.4365, 0.4370, 0.4375, 0.4380, 0.4385, 0.4390, 0.4395, 0.44,
                0.4405, 0.4410, 0.4415, 0.4420, 0.4425, 0.4430, 0.4435, 0.4440, 0.4445, 0.4450, 0.4455, 0.4460, 0.4465, 0.4470, 0.4475,
                0.4480, 0.4485, 0.4490, 0.4495, 0.45, 0.4505, 0.4510, 0.4515, 0.4520, 0.4525, 0.4530, 0.4535, 0.4540, 0.4545, 0.4550,
                0.4555, 0.4560, 0.4565, 0.4570, 0.4575, 0.4580, 0.4585, 0.4590, 0.4595, 0.46]"""

temperatures = [0.5, 0.502, 0.504, 0.506, 0.508, 0.51, 0.512, 0.514, 0.516, 0.518, 0.52, 
                0.522, 0.524, 0.526, 0.528, 0.53, 0.532, 0.534, 0.536, 0.538, 0.54, 0.542, 
                0.544, 0.546, 0.548, 0.55, 0.552, 0.554, 0.556, 0.558, 0.56, 0.562, 0.564, 
                0.566, 0.568, 0.57, 0.572, 0.574, 0.576, 0.578, 0.58, 
                0.582, 0.584, 0.586, 0.588, 0.59, 0.592, 0.594, 0.596, 0.598, 0.6]

temperatures = [ 
     0.4, 0.45, 0.5] 

#L = [8, 16, 32, 64]
#L = [8, 12, 16, 20, 24, 32, 40, 64 ]
L=[8]

K=1
lines = []

def calculate_std(data): 
    try:
        std_deviation = np.std(data)
        return std_deviation
    except Exception as e:
        print(f"Error: {e}")
        return None

def calculate_mean(data): 
    mean = sum(data) / len(data)
    return mean

i = 0 

for l in L : 

    N= l * l

    pseudo_magn = []
    mm_2_values = []
    mm_4_values = [] 
    cumulant = []   

    for temp in temperatures: 

        #file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output_K={K}/Output_L={l}/T_{temp}" + '/trsb_magnetization.txt'
        #file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output_prova/Output_L={l}/T_{temp}" + '/trsb_magnetization.txt'
        file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Relative_fluctuations/K={K}/Output_L={l}/T_{temp}" + '/trsb_magnetization.txt'


        with open(file_path, 'r') as file:
            numbers = [float(line.strip()) for line in file.readlines()]
            m2= [x**2 for x in numbers ]
            m4= [x**4 for x in numbers ]
            avg_m2 = calculate_mean (m2)
            avg_m4 = calculate_mean (m4)
            mm = calculate_mean(numbers)
            U = avg_m4 / (3 * avg_m2**2)
            cumulant.append (np.abs(U))
            pseudo_magn.append(np.abs(mm))
    


    #Cumulant Binder 
    plt.plot(temperatures, cumulant, linestyle='-', label = f'L={l}')
    cum_array = np.array (cumulant) 
    lines.append (cum_array)



plt.xlabel('T(K)')
plt.ylabel('U')
plt.title('Relative density fluctuations ')
plt.legend()
plt.grid(True)





#INTERSECTION 
intersections = []

lin1 = lines [0] 
lin2 = lines [1]

Line1 = LineString (np.column_stack((temperatures, lin1)))
Line2 = LineString (np.column_stack((temperatures, lin2)))
intersection = Line1.intersection (Line2)
print (intersection)
intersections.append(intersection)

plt.plot(*intersection.xy, 'ro')

plt.show()