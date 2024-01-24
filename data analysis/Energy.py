import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 

temperatures = [ 0.009, 0.01, 0.02, 0.03, 0.05, 0.07, 0.08, 0.09, 0.1, 0.11, 0.13, 0.14, 0.15, 
                0.17, 0.2, 0.22, 0.25, 0.27, 0.3, 0.32, 0.35, 0.37, 0.4, 0.45, 0.5, 0.52, 0.55, 
                0.57, 0.6, 0.62, 0.65, 0.67, 0.7, 0.72, 0.75, 0.77, 0.8, 0.82, 0.85, 0.87, 0.9, 
                0.92, 0.94, 0.95, 0.96, 0.98, 1.0, 1.02, 1.04, 1.05, 1.06, 1.07, 1.09, 1.1, 1.13, 1.15, 1.17, 1.2, 1.22, 1.25, 
                1.27, 1.3, 1.4, 1.45, 1.5, 1.7, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5]  

L= 20
K=11

N = L * L

# Create a dictionary to store the arrays
mean_argument_values = []
std_dev = []

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



#PLOT OF THE ENERGY 
i = 0
for temp in temperatures:
        

    file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output_K={K}/Output_L={L}/T_{temp}" + '/Energy.txt'

    with open(file_path, 'r') as file:
        numbers = [float(line.strip()) for line in file.readlines()]
        mm = calculate_mean(numbers)
        en_std = calculate_std(numbers)
        mean_argument_values.append(mm) 
        std_dev.append(en_std)
    
mean_val = np.array (mean_argument_values)
std_val = np.array (std_dev)

# Plot Energy vs. Temperature
plt.plot(temperatures, mean_val, linestyle='-')  #marker='o',
plt.fill_between(temperatures, mean_val - std_val, mean_val + std_val, alpha=0.6, facecolor='lightblue', linewidth=4, label='Error Bounds ')

plt.xlabel('Temperature (K)')
plt.ylabel('Energy')
plt.title(f'Mean Energy  ($J_1 =J_2, K=2$), L={L} ')
plt.grid(True)

plt.show()