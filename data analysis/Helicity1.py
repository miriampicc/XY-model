import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 

temperatures = [0.5, 0.502, 0.504, 0.506, 0.508, 0.51, 0.512, 0.514, 0.516, 0.518, 0.52, 
                0.522, 0.524, 0.526, 0.528, 0.53, 0.532, 0.534, 0.536, 0.538, 0.54, 0.542, 
                0.544, 0.546, 0.548, 0.55, 0.552, 0.554, 0.556, 0.558, 0.56, 0.562, 0.564, 
                0.566, 0.568, 0.57, 0.572, 0.574, 0.576, 0.578, 0.58, 
                0.582, 0.584, 0.586, 0.588, 0.59, 0.592, 0.594, 0.596, 0.598, 0.6]

L=40
K=5

N = L * L
# Create a dictionary to store the arrays
mean_argument_values = []

def calculate_std(data): 
    try:
        std_deviation = np.std(data)
        return std_deviation
    except Exception as e:
        # Handle any exceptions that may occur during the calculation
        print(f"Error: {e}")
        return None

def calculate_mean(data): 
    mean = sum(data) / len(data)
    return mean


# Superfluid Stiffness
Jp1 = []
Js1 = []
mean_Jd_values1 = [] 
sin_Ic1= []

i = 0 

for temp in temperatures:
    Jd1 = []
    Ic1 = []

    file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Relative_fluctuations/K={K}/Output_L={L}/T_{temp}" + '/Helicity_modulus1.txt'
    #file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Relative_fluctuations/K={K}/Output_L={L}/T_{temp}" + '/Helicity_modulus1.txt'

    with open(file_path, "r") as file:
        for line in file:
            columns = line.strip().split()
            if len(columns) == 2:
                Jd1.append(float(columns[0]))
                Ic1.append(float(columns[1]))

    mean_sin = calculate_std(Ic1) ** 2
    sin_Ic1.append(mean_sin)
    sin = mean_sin * N / temperatures[i]

    Jp1.append(sin)

    mm = calculate_mean(Jd1) 
    cos_Jd = mm 
    mean_Jd_values1.append(cos_Jd)

    Js_new = cos_Jd - sin
    Js1.append(Js_new)
    
    i += 1 


Jd_matrix1 = np.column_stack((temperatures, mean_Jd_values1))
Jp_matrix1 = np.column_stack((temperatures, Jp1))
Js_matrix1 = np.column_stack((temperatures, Js1))

y = [0] * len(temperatures)

for temp in range(len(temperatures)): 
    y[temp] = 2 * temperatures[temp] / np.pi

plt.plot(Jd_matrix1[:, 0], Jd_matrix1[:, 1], linestyle='-', label=r'$J_d$')
plt.plot(Jp_matrix1[:, 0], Jp_matrix1[:, 1],  linestyle='-', label = r'$J_p$')
plt.plot(Js_matrix1[:, 0], Js_matrix1[:, 1], linestyle='-', label = r'$J_s$')
plt.plot(temperatures , y, linestyle='--', label = r'$\frac{2T}{\pi}$')
plt.legend()
plt.xlabel('Temperature (K)')
plt.ylabel('superfluid stiffness')
plt.title('Superfulitity Stiffness Lattice 1, ($J_1 =2 J_2=K$), n$_{steps}$= 100000')
plt.grid(True)

plt.show()