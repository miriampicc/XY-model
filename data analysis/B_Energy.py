import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 
import argparse

parser = argparse.ArgumentParser(description='Description of the script')

# Add arguments
parser.add_argument('--L', nargs='+', type=int, help='L')
parser.add_argument('--K', type=float, help='K')
parser.add_argument('--e', type=float, help='e')
parser.add_argument('--b_high', type=float, help='beta high')
parser.add_argument('--b_low', type=float, help='beta low')
parser.add_argument('--rank', type=int, help='rank')

# Parse the command-line arguments
args = parser.parse_args()

# Access the argument values
L = args.L
K = args.K
e = args.e
beta_high = args.b_high
beta_low = args.b_low
rank = args.rank

# Now you can use these values in your script
print("L=", L)
print("K=", K)
print("e=", e)
print("beta high=", beta_high)
print("beta low=", beta_low)
print("rank=", rank)



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
    N = l * l

    for n in range(rank) : 

        t = T_high - n * delta
        #print (t)
        
        temperatures.append(t)

        #file_path = f"/Users/mirimi/Desktop/OUTPUT_cluster/e_{e}/L{l}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}/beta_{n}" + '/Energy.txt'
        file_path = f"/home/x_mirpi/Output_TBG/K_{K}_first/e_{e}/L{l}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}/beta_{n}" + '/B_Energy.txt'
        #print(n)

        with open(file_path, 'r') as file:
            energies = [float(line.strip()) for line in file.readlines()]
            #print(len(energies))
            en_avg = calculate_mean(energies)
            en_std = calculate_std(energies)
            mean_energy_values.append(en_avg) 
            std_dev.append(en_std)

    
    energy = np.array(mean_energy_values)
    energy = energy / N
    temp = np.array(temperatures) 
    std_val = np.array(std_dev)

    # Plot Energy vs. Temperature
    plt.plot(temp, energy, linestyle='-', label=f'L={l}')  #marker='o',
    plt.fill_between(temp, energy - std_val, energy + std_val, alpha=0.3, linewidth=4)

#plt.ylim(top=0)
plt.xlabel('Temperature (K)')
plt.ylabel(' Energy ')
plt.title(f'B Energy e={e}, K={K}')
plt.legend()
plt.grid(True)
plt.savefig(f'B_Energy_e={e}_K={K}.jpg')

plt.show()