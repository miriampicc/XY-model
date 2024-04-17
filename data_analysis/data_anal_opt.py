import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 

temperatures = [ 0.009, 0.01, 0.02, 0.03, 0.05, 0.07, 0.08, 0.09, 0.1, 0.11, 0.13, 0.14, 0.15, 
                0.17, 0.2, 0.22, 0.25, 0.27, 0.3, 0.32, 0.35, 0.37, 0.4, 0.45, 0.5, 0.52, 0.55, 
                0.57, 0.6, 0.62, 0.65, 0.67, 0.7, 0.72, 0.75, 0.77, 0.8, 0.82, 0.85, 0.87, 0.9, 
                0.92, 0.94, 0.95, 0.96, 0.98, 1.0, 1.02, 1.04, 1.05, 1.06, 1.07, 1.09, 1.1, 1.13, 
                1.15, 1.17, 1.2, 1.22, 1.25, 1.27, 1.3, 1.4, 1.45, 1.5, 1.7, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5]  
values = []

mean_vertex_values1 = []
mean_vertex_values2 = []
mean_antivertex_values1 = []
mean_antivertex_values2 = []

N = 20 * 20

def calculate_std (data) : 
    try:

        std_deviation = np.std(data)
        return std_deviation
    
    except Exception as e:
        # Handle any exceptions that may occur during the calculation
        print(f"Error: {e}")
        return None


def calculate_mean (data) : 

    mean = np.mean(data)
    return mean

#Vortices and Antivortices 

for temp in temperatures  :
    layer1_values = []
    layer2_values = []

    file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output_L=20/T_{temp}" + '/Vortices.txt'
    
    with open(file_path, 'r') as file:

        for line in file: 

            values = [float(val) for val in line.strip().split()]
            layer1_values.append(values[0])
            layer2_values.append(values[1])

        v1 = calculate_mean (layer1_values)
        v2 = calculate_mean (layer2_values)

        mean_vertex_values1.append(v1)
        mean_vertex_values1.append(v2)
            
vertex_matrix1 = np.column_stack((temperatures, mean_vertex_values1))
vertex_matrix2 = np.column_stack((temperatures, mean_vertex_values2))

for temp in temperatures  :
    layer1_values = []
    layer2_values = []

    file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output_L=20/T_{temp}" + '/Antivortices.txt'

    with open(file_path, 'r') as file:

        for line in file: 

            values = [float(val) for val in line.strip().split()]
            layer1_values.append(values[0])
            layer2_values.append(values[1])

        v1 = calculate_mean (layer1_values)
        v2 = calculate_mean (layer2_values)

        mean_antivertex_values1.append(v1)
        mean_antivertex_values1.append(v2)
            
antivertex_matrix1 = np.column_stack((temperatures, mean_antivertex_values1))
antivertex_matrix2 = np.column_stack((temperatures, mean_antivertex_values2))



# Create a new figure for the first set of plots
fig1, ax1 = plt.subplots()
ax1.plot(vertex_matrix1[:, 0], vertex_matrix1[:, 1], marker='o', linestyle='-', label='Vortex Layer 1')
ax1.plot(antivertex_matrix1[:, 0], antivertex_matrix1[:, 1], marker='o', linestyle='-', label='Antivortex Layer 1')
ax1.legend()
ax1.set_xlabel('Temperature (K)')
ax1.set_ylabel('Vortex, Antivortex')
ax1.set_title('Vortices and Antivortices vs. Temperatures (Layer 1, L=16, J=1)')
ax1.grid(True)

# Create a new figure for the second set of plots
fig2, ax2 = plt.subplots()
ax2.plot(vertex_matrix2[:, 0], vertex_matrix2[:, 1], marker='o', linestyle='-', label='Vortex Layer 2')
ax2.plot(antivertex_matrix2[:, 0], antivertex_matrix2[:, 1], marker='o', linestyle='-', label='Antivortex Layer 2')
ax2.legend()
ax2.set_xlabel('Temperature (K)')
ax2.set_ylabel('Vortex, Antivortex')
ax2.set_title('Vortices and Antivortices vs. Temperatures (Layer 2, L=16, J=1)')
ax2.grid(True)

# Display the plots
plt.show()

