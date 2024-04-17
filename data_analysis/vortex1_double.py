import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 


##temperatures = [ 0.009 , 0.01 , 0.02 , 0.03 , 0.05 , 0.07 , 0.08 , 0.09 , 0.1 , 0.11 , 0.13 , 0.14 , 0.15, 0.2 , 0.25 , 0.3 , 0.35 , 
  #              0.4 , 0.45 , 0.5 , 0.55 , 0.6 , 0.65 , 0.7 , 0.75 , 0.8 , 0.85 , 0.9 , 0.92 , 0.94 , 0.95 , 0.96 , 0.98 , 1.0 , 1.02 , 1.04 , 1.05 , 1.07 ,
   #               1.09 , 1.1 , 1.12 , 1.14 , 1.15 , 1.17 , 1.18 , 1.2 , 1.23 , 1.25 , 1.27 , 1.29 , 
    #            1.3 , 1.33 , 1.35 , 1.37 , 1.4 , 1.45 , 1.5 , 1.6 , 1.7 , 1.8 , 1.9 , 2.0 , 2.1 , 2.2 , 2.3 , 2.4 , 2.5 , 
     #           3.0 , 3.5 , 4.0 , 4.5 , 5.0  ]  


temperatures = [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,
                0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.27, 1.3, 1.4, 1.45, 1.5 ]

N = 16 * 16

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


mean_vertex_values = []
mean_antivertex_values = []

i = 0 

for temp in temperatures:
    file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output/T_{temp}" + '/Vortices2.txt'
    
    with open(file_path, 'r') as file:
        numbers = [float(line.strip()) for line in file.readlines()]
        mm = calculate_mean(numbers)
        mean_vertex_values.append(mm)
    
    print("vortices ", i)
    i += 1
            
vertex_matrix = np.column_stack((temperatures, mean_vertex_values))

i = 0

mean_antivertex_values = []

for temp in temperatures:
    file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output/T_{temp}" + '/Antivortices2.txt'

    with open(file_path, 'r') as file:
        numbers = [float(line.strip()) for line in file.readlines()]
        mm = calculate_mean(numbers)
        mean_antivertex_values.append(-mm)

    print("antivortices ", i)
    i += 1
            
antivertex_matrix = np.column_stack((temperatures, mean_antivertex_values))


# Plot Vortex and Antivortex vs. Temperature
plt.plot(vertex_matrix[:, 0], vertex_matrix[:, 1], marker='o', linestyle='-', label='Vortex')
plt.plot(antivertex_matrix[:, 0], antivertex_matrix[:, 1], marker='o', linestyle='-', label='Antivortex')
plt.legend()
plt.xlabel('Temperature (K)')
plt.ylabel('Vortex, Antivortex')
plt.title('Vortices and Antivortices (Lattice 2) ($J_1 = J_2 K = 1.5$), $n_{steps}= 10000$')
plt.grid(True)

plt.show()