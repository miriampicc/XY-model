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
parser.add_argument('--beta', type=int, help='beta')
parser.add_argument('--num_bins', type=int, help='num_bins')


# Parse the command-line arguments
args = parser.parse_args()

# Access the argument values
L = args.L
K = args.K
e = args.e
beta_high = args.b_high
beta_low = args.b_low
rank = args.rank
beta = args.beta
beta = args.beta
num_bins = args.num_bins

# Now you can use these values in your script
print("L=", L)
print("K=", K)
print("e=", e)
print("beta high=", beta_high)
print("beta low=", beta_low)
print("rank=", rank)
print("beta=", beta)
print("number of bins=", num_bins)


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

colors = [
    "#00008B",  # Dark Blue
    "#0000CD",  # Medium Blue
    "#4169E1",  # Royal Blue
    "#1E90FF",  # Dodger Blue
    "#00BFFF",  # Deep Sky Blue
    "#87CEEB",  # Sky Blue
    "#87CEFA",  # Light Sky Blue
    "#EEE8AA",  # Pale Goldenrod
    "#F0E68C",  # Khaki
    "#FFD700",  # Gold
    "#F08080",  # Light Coral
    "#FF6347"   # Tomato
]

#let us obtain the temperatures 
T = 1/beta

i = 0

for l in L:

    N = l * l

    file_path = f"/home/x_mirpi/Output_TBG/K_{K}_tdf/e_{e}/L{l}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}/beta_{beta}" + '/Energy.txt'

    with open(file_path, 'r') as file:
        energies = [float(line.strip()) for line in file.readlines()]
            #print(len(energies))
    steps = len(energies)
    energy = np.array(energies)
    energy = energy / N

    log_bins = np.logspace(0, np.log10(steps), num_bins, base=10.0) 

    mean_values = []
    std_dev_values = []
    half_widths = []
    inside_bin = []
    x_mc_val = []

    i = 0 
    j = 0


    for j in range(len(log_bins)-1) : 
        inside_bin = []

        for i in range(len(energy)):
            if log_bins[j] <= i < log_bins[j + 1]:
                inside_bin.append(energy[i])

        mean_val = np.mean(inside_bin)
        std_val = np.std(inside_bin)
        x_middle = (log_bins[j + 1] + log_bins[j]) / 2.0

        mean_values.append(mean_val)
        std_dev_values.append(std_val)
        half_widths.append(x_middle)

    print (half_widths)


#plt.ylim(top=0)
plt.scatter(half_widths, mean_values, label='Mean Energy ', marker='o')
plt.errorbar(half_widths, mean_values, yerr=std_dev_values,linestyle='None', color='black', capsize=3, label='Error bars')
plt.xscale('log')  # Set the x-axis to be logarithmic
plt.xlabel('Monte Carlo Steps (log scale)')
plt.ylabel('Energy')
plt.legend()
plt.title(f'Energy Thermalization, T={T}, MC ={steps}, L={l} ')
plt.grid()
plt.savefig(f'T_Energy_e={e}_K={K}_L{l}.jpg')
plt.show()
