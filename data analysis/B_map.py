import numpy as np 
import matplotlib as plt 

import argparse

parser = argparse.ArgumentParser(description='Description of the script')

# Add arguments
parser.add_argument('--L', type=int, help='L')
parser.add_argument('--K', type=float, help='K')
parser.add_argument('--e', type=float, help='e')
parser.add_argument('--beta', type=int, help='beta')
parser.add_argument('--b_high', type=float, help='beta high')
parser.add_argument('--b_low', type=float, help='beta low')
parser.add_argument('--rank', type=int, help='rank')
parser.add_argument('--N_steps', type=float, help='N_steps' )

# Parse the command-line arguments
args = parser.parse_args()

# Access the argument values
L = args.L
K = args.K
e = args.e
beta = args.beta
beta_high = args.b_high
beta_low = args.b_low
rank = args.rank
N_steps = args.N_steps

# Now you can use these values in your script
print("L=", L)
print("K=", K)
print("e=", e)
print("beta =", beta)
print("beta high=", beta_high)
print("beta low=", beta_low)
print("rank=", rank)
print("Number of steps =", N_steps)


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


N = L * L

file_path = f"/home/x_mirpi/Output_TBG/K_{K}_sB/e_{e}/L{L}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}/beta_{beta}" + '/B_Plaquette.txt'

with open(file_path, 'r') as file:
    all_values = [float(line.strip()) for line in file.readlines()]

lattices = []

# Divide data into lattices of size N
lattices = [all_values[i:i + N] for i in range(0, len(all_values), N)]

    
mean_energy = []

for site in range(N):

    Plaquette_en = []

    for n_array in range(N_steps): 

        val = lattices [n_array][site]
        Plaquette_en.append(val)

    plaqette = np.array(Plaquette_en)

    energy_avg = calculate_mean(plaqette)
    mean_energy.append(energy_avg)


energy = np.array(mean_energy)
energy = energy / N

# Reshape energy into a 2D array of size LxL
energy_matrix = energy.reshape(L, L)

plt.figure(figsize=(8, 6))
plt.imshow(energy_matrix, cmap='viridis', interpolation='nearest', origin='lower')
plt.colorbar(label='Mean Energy')
plt.title('Mean Energy Map')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.grid(False)  # Disable grid lines
plt.savefig(f'Color_map_L={L}_beta={beta}.jpg')
plt.show()
