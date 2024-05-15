import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 

L = [8, 16, 24, 32, 48, 64]
K = 5.0
e = 0.1
beta_high = 1.786   #T=0.56  1.786  #1.754
beta_low = 1.695   #T=0.59
rank = 64


# Create a dictionary to store the arrays


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


#let us obtain the temperatures 
temperatures = []
delta = (1/beta_low - 1/beta_high)/(rank)
T_high = 1/beta_low
T_low = 1/beta_high
print(delta)

for l in L:
    mean_energy_values = []
    std_dev = []
    temperatures = []

    for n in range(rank) : 

        t = T_high - n * delta
        #print (t)
        
        temperatures.append(t)

        file_path = f"/home/x_mirpi/Output_TBG/K_{K}_first/e_{e}/L{l}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}/beta_{n}" + '/Energy.txt'
        #print(n)

        with open(file_path, 'r') as file:
            energies = [float(line.strip()) for line in file.readlines()]
            #print(len(energies))
            en_avg = calculate_mean(energies)
            en_std = calculate_std(energies)
            mean_energy_values.append(en_avg) 
            std_dev.append(en_std)

    
    energy = np.array(mean_energy_values)
    temp = np.array(temperatures) 
    std_val = np.array(std_dev)

    # Plot Energy vs. Temperature
    plt.plot(temp, energy, linestyle='-', label=f'L={l}')  #marker='o',
    plt.fill_between(temp, energy - std_val, energy + std_val, alpha=0.6, linewidth=4)

#plt.ylim(top=0)
plt.xlabel('Temperature (K)')
plt.ylabel(' Energy ')
plt.title(f'Total Energy')
plt.legend()
plt.grid(True)

plt.show()