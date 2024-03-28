import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 

numbers = ['1', '2']
arguments = ['/Magnetization1.txt', '/Magnetization2.txt'] 
title = ['Lattice 1', 'Lattice 2']

temperatures = [ 0.009, 0.01, 0.02, 0.03, 0.05, 0.07, 0.08, 0.09, 0.1, 0.11, 0.13, 0.14, 0.15, 
                0.17, 0.2, 0.22, 0.25, 0.27, 0.3, 0.32, 0.35, 0.37, 0.4, 0.45, 0.5, 0.52, 0.55, 
                0.57, 0.6, 0.62, 0.65, 0.67, 0.7, 0.72, 0.75, 0.77, 0.8, 0.82, 0.85, 0.87, 0.9, 
                0.92, 0.94, 0.95, 0.96, 0.98, 1.0, 1.02, 1.04, 1.05, 1.06, 1.07, 1.09, 1.1, 1.13, 1.15, 1.17, 1.2, 1.22, 1.25, 
                1.27, 1.3, 1.4, 1.45, 1.5, 1.7, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5] 

L=[8, 12, 16, 20, 24, 32, 40]
rainbow_colors = ['red', 'orange', 'yellow', 'green', 'blue', 'indigo', 'violet']

K=1


def calculate_std(data): 
    try:
        std_deviation = np.std(data)
        std_deviation = std_deviation / (np.sqrt(len(data)-1))
        return std_deviation
    except Exception as e:
        # Handle any exceptions that may occur during the calculation
        print(f"Error: {e}")
        return None

def calculate_mean(data): 
    mean = sum(data) / len(data)
    return mean


i =0

for l in L: 
    N = l * l

    mean_argument_values = {0: [], 1: []}
    std_values = {0: [], 1: []}
 
    for temp in temperatures:
        layer1_values = []
        layer2_values = []
        mm1=0
        mm2=0
        std1=0
        std2=0

        #file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Uppsala/J={j}/Output_L={L}/T_{temp}" + '/Single_Magnetization.txt'
        #file_path = f"/Users/mirimi/Desktop//XY-model-single/J={J}/Output_L={l}/T_{temp}" + '/Single_Magnetization.txt'
        #file_path = f"/Users/mirimi/Desktop//XY-model-single/Prova2/Output_L={L}/T_{temp}" + '/Single_Magnetization.txt'
        file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Relative_fluctuations/K={K}/Output_L={l}/T_{temp}" + '/Single_Magnetization.txt'




        with open(file_path, 'r') as file:

            for line in file:
                    # Assuming the file has two columns separated by whitespace
                    values = [float(val) for val in line.strip().split()]
                    layer1_values.append(values[0])
                    layer2_values.append(values[1])

            print(len(layer1_values))

            mm1 = calculate_mean(layer1_values)
            mm2 = calculate_mean(layer2_values)
            std1 = calculate_std (layer1_values)
            std2 = calculate_std (layer2_values)

            mean_argument_values[0].append(mm1)
            mean_argument_values[1].append(mm2)
            std_values[0].append(std1)
            std_values[1].append(std2)        
    
    magn1 = np.array(mean_argument_values[0])
    magn2 = np.array (mean_argument_values[1])
    m_std1 = np.array(std_values[0])
    m_std2 = np.array(std_values[1])

    plt.plot(temperatures, mean_argument_values[0])
    plt.fill_between(temperatures, magn1 - m_std1, magn1 + m_std1, alpha=0.3, linewidth=4)

    plt.plot(temperatures, mean_argument_values[1], label=f'L= {l}', color=rainbow_colors[i % len(rainbow_colors)])
    plt.fill_between(temperatures, magn2 - m_std2, magn2 + m_std2, alpha=0.2, linewidth=3, color=rainbow_colors[i % len(rainbow_colors)])

    i +=1


plt.xlabel('Temperature')
plt.ylabel('Magnetization Values')
plt.title(f'Relative desity fluctations, K={K}')
plt.legend()
plt.grid(True)

plt.show()