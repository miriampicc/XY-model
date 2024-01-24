import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 

temperatures = [ 0.009, 0.01, 0.02, 0.03, 0.05, 0.07, 0.08, 0.09, 0.1, 0.11, 0.13, 0.14, 0.15, 
                0.17, 0.2, 0.22, 0.25, 0.27, 0.3, 0.32, 0.35, 0.37, 0.4, 0.45, 0.5, 0.52, 0.55, 
                0.57, 0.6, 0.62, 0.65, 0.67, 0.7, 0.72, 0.75, 0.77, 0.8, 0.82, 0.85, 0.87, 0.9, 
                0.92, 0.94, 0.95, 0.96, 0.98, 1.0, 1.02, 1.04, 1.05, 1.06, 1.07, 1.09, 1.1, 1.13, 1.15, 1.17, 1.2, 1.22, 1.25, 
                1.27, 1.3, 1.4, 1.45, 1.5, 1.7, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5]  

L=8
K=0.5

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

    file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output_K={K}/Output_L={L}/T_{temp}" + '/Helicity_modulus1.txt'
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