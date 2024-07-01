import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
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

temperatures = []
delta = (1/beta_low - 1/beta_high)/(rank)
T_high = 1/beta_low
T_low = 1/beta_high
print(delta)

delta_beta = (beta_high - beta_low)/(rank)

colors = [
    "#00008B",  # Dark Blue
    "#0000CD",  # Medium Blue
    "#1E90FF",  # Dodger Blue
    "#87CEEB",  # Sky Blue
    "#EEE8AA",  # Pale Goldenrod
    "#F0E68C",  # Khaki
    "#FFD700",  # Gold
    "#F08080",  # Light Coral
    "#FF6347",  # Tomato
    "#32CD32"   # Lime Green
]

# PLOT OF THE SPECIFIC HEAT
N_dataset = 1000
block_size = 1000
i = 0
for l in L:
    N = l * l

    specific_heat = []
    betas = []
    temperatures = []
    val = []
    error = []

    for n in range(rank):
        t = T_high - n * delta
        bb = beta_low + delta_beta * n

        temperatures.append(t)
        betas.append(bb)

        file_path = f"/home/x_mirpi/Output_TBG/K_{K}/e_{e}/L{l}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}/beta_{n}" + '/Energy.txt'

        with open(file_path, 'r') as file:
            numbers = []
            for line_number, line in enumerate(file, start=1):
                stripped_line = line.strip()
                try:
                    number = float(stripped_line)
                    numbers.append(number)
                except ValueError:
                    print(f"Warning: Could not convert line {line_number} to float: {stripped_line}")

        blocks = [numbers[j:j + block_size] for j in range(0, len(numbers), block_size)]

        specific_heat = []

        for n_boots in range(N_dataset):
            new_dataset = []

            for b in blocks:
                rand_index = random.randint(1, len(blocks)-1)
                new_dataset.append(blocks[rand_index])

            energy_dataset = np.array([])

            for block in new_dataset:
                energy_dataset = np.append(energy_dataset, block)  # Append each block to the dataset

            mm = calculate_mean(energy_dataset)
            en_var = (np.var(energy_dataset))
            cc = (((1/t)**2) * en_var )/ N

            specific_heat.append(cc)

        sh_val = np.array(specific_heat)
        sh = calculate_mean(sh_val)
        sh_err = np.std(sh_val)

        val.append(sh)
        error.append(sh_err)

        beta_array = np.array(betas)

    # Plot Cv vs. Temperature
    plt.plot(betas, val, linestyle='-', label=f'L={l}', color=colors[i])
    plt.fill_between(betas, val - error, val + error, color=colors[i], alpha=0.3)

    i = i + 1

plt.xlabel(r'$\beta$')
plt.ylabel('$C_V$')
plt.title(f'Specific Heat $K = {K}$, $e={e}$ ')
plt.legend()
plt.grid(True)
plt.savefig(f'Specific_heat_b_e={e}_K={K}.jpg')

plt.show()
