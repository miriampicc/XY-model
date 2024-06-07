import numpy as np
from statsmodels.tsa.stattools import acf
import matplotlib.pyplot as plt
import argparse
import os

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

# Print arguments to confirm
print("L=", L)
print("K=", K)
print("e=", e)
print("beta high=", beta_high)
print("beta low=", beta_low)
print("rank=", rank)

# Delta calculation
delta = (1/beta_low - 1/beta_high) / rank
T_high = 1/beta_low
T_low = 1/beta_high
print("Delta:", delta)

# Prepare lists for storing results
temperatures = []
correlation_times1 = []

for l in L:
    autocorrelations_layer1 = []

    for n in range(rank):
        t = T_high - n * delta
        temperatures.append(t)

        file_path = f"/home/x_mirpi/Output_TBG/K_{K}/e_{e}/L{l}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}/beta_{n}/Single_Magnetization.txt"

        m1 = []
        m2 = []

        # Read data from file
        if os.path.exists(file_path):
            with open(file_path, "r") as file:
                for line in file:
                    columns = line.strip().split()
                    if len(columns) == 2:
                        m1.append(float(columns[0]))
                        m2.append(float(columns[1]))
        else:
            print(f"File not found: {file_path}")
            continue

        # Calculate autocorrelations
        autocorrelation_1 = acf(m1, nlags=50)
        autocorrelations_layer1.append(autocorrelation_1)

        # Calculate \tau by finding the lag where autocorrelation drops below a threshold
        threshold = 1 / 2.71828  # 1/e
        tau = None
        for lag, acf_value in enumerate(autocorrelation_1):
            if acf_value < threshold:
                tau = lag
                break
        correlation_times1.append(tau)

    lags1 = np.arange(len(autocorrelations_layer1[0]))

    # Plot all autocorrelation functions on the same plot
    plt.figure(figsize=(10, 6))
    for i, autocorr in enumerate(autocorrelations_layer1):
        label = f'T={temperatures[i]:.2f}'
        plt.plot(lags1, autocorr, linestyle='-', label=label) #marker='o',

plt.title('Autocorrelation Functions for Different Temperatures')
plt.xlabel('Lag')
plt.ylabel('Autocorrelation')
plt.grid(True)
plt.savefig(f'Autocorr_lag_K{K}_e{e}.jpg')
plt.legend()
plt.show()


# Plot the correlation times as a function of temperature
plt.figure(figsize=(10, 6))
plt.plot(temperatures, correlation_times1, linestyle='-') #marker='o', 
plt.title('Correlation Time (\u03C4) as a Function of Temperature')
plt.xlabel('Temperature')
plt.ylabel('Correlation Time (\u03C4)')
plt.savefig(f'Autocorr_temp_K{K}_e{e}.jpg')
plt.grid(True)
plt.show()
