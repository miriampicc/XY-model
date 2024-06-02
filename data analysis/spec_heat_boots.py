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


def calculate_mean(data): 
    mean = sum(data) / len(data)
    return mean

def bootstrap(data, statistic, n_bootstrap=1000):
    n = len(data)
    bootstrap_estimates = np.empty(n_bootstrap)
    
    for i in range(n_bootstrap):
        bootstrap_sample = np.random.choice(data, size=n, replace=True)
        bootstrap_estimates[i] = statistic(bootstrap_sample)
    
    return bootstrap_estimates

def specific_heat_statistic(data, t, N):
    en_var = (np.std(data))**2
    cc = (((1/t)**2) * en_var )/ N
    return cc

temperatures = []
delta = (1/beta_low - 1/beta_high)/(rank)
T_high = 1/beta_low
T_low = 1/beta_high
print(delta)

delta_beta = (beta_high - beta_low)/(rank)



#PLOT OF THE SPECIFIC HEAT 
i = 0
for l in L:
    N = l * l

    specific_heat = []
    betas = []
    temperatures = []
    specific_heat_errors = []


    for n in range(rank):

        t = T_high - n * delta
        bb = beta_low + delta_beta * n 
        
        temperatures.append(t)
        betas.append(bb)
        

        file_path = f"/home/x_mirpi/Output_TBG/K_{K}/e_{e}/L{l}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}/beta_{n}" + '/Energy.txt'

        with open(file_path, 'r') as file:
            numbers = [float(line.strip()) for line in file.readlines()]

            bootstrap_estimates = bootstrap(numbers, lambda x: specific_heat_statistic(x, t, N))
            mean_specific_heat = np.mean(bootstrap_estimates)
            std_error_specific_heat = np.std(bootstrap_estimates)

            specific_heat.append(mean_specific_heat)
            specific_heat_errors.append(std_error_specific_heat)

        
    sh_val = np.array(specific_heat)
    sh_err = np.array(specific_heat_errors)
    beta_array = np.array(betas)

    # Plot Energy vs. Temperature with error bars
    plt.errorbar(betas, sh_val, yerr=sh_err, linestyle='-', label=f'L={l}')
   

plt.xlabel(r'$\beta$')
plt.ylabel('$C_V$')
plt.title(f'Specific Heat $K = {K}$, $e={e}$ ')
plt.legend()
plt.grid(True)
plt.savefig(f'Specific_heat_e={e}_K={K}.jpg')


plt.show()