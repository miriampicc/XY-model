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
parser.add_argument('--a', type=float, help='rank')


# Parse the command-line arguments
args = parser.parse_args()

# Access the argument values
L = args.L
K = args.K
e = args.e
beta_high = args.b_high
beta_low = args.b_low
rank = args.rank
a = args.a

# Now you can use these values in your script
print("L=", L)
print("K=", K)
print("e=", e)
print("beta high=", beta_high)
print("beta low=", beta_low)
print("rank=", rank)
print("a=", a)

def calculate_mean(data):
    mean = sum(data) / len(data)
    return mean

temperatures = []
delta = (1 / beta_low - 1 / beta_high) / rank
T_high = 1 / beta_low
T_low = 1 / beta_high
print(delta)

delta_beta = (beta_high - beta_low) / rank

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

# PLOT OF THE SPECIFIC HEAT
i = 0
for l in L:
    N = l * l

    specific_heat = []
    betas = []
    temperatures = []

    for n in range(rank):
        t = T_high - n * delta
        bb = beta_low + delta_beta * n

        temperatures.append(t)
        betas.append(bb)

        #file_path = f"/home/x_mirpi/Output_TBG/K_{K}_tdf2/e_{e}/L{l}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}/beta_{n}/Energy.txt"
        file_path = f"/home/x_mirpi/Output_TBG/K_{K}_tdf2/e_{e}/L{l}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}_a{a}/beta_{n}" + '/Energy.txt'


        try:
            with open(file_path, 'r') as file:
                numbers = []
                for line in file:
                    line = line.strip()
                    if line:
                        try:
                            numbers.append(float(line))
                        except ValueError:
                            print(f"Warning: Could not convert line to float: '{line}'")
                            continue
                if numbers:
                    mm = calculate_mean(numbers)
                    en_var = (np.std(numbers)) ** 2
                    cc = (((1 / t) ** 2) * en_var) / N
                    specific_heat.append(cc)
                else:
                    specific_heat.append(np.nan)  # Append NaN if no valid data
        except FileNotFoundError:
            print(f"File not found: {file_path}")
            specific_heat.append(np.nan)  # Append NaN if file not found

    sh_val = np.array(specific_heat)
    beta_array = np.array(betas)

    # Ensure both arrays have the same length
    if len(beta_array) == len(sh_val):
        # Plot Energy vs. Temperature
        plt.plot(beta_array, sh_val, linestyle='-', label=f'L={l}', color=colors[i])
        i += 1
    else:
        print(f"Length mismatch: betas ({len(beta_array)}) and specific_heat ({len(sh_val)})")

plt.xlabel(r'$\beta$')
plt.ylabel('$C_V$')
plt.title(f'Specific Heat $K = {K}$, $a={a}$ ')
plt.legend()
plt.grid(True)
plt.savefig(f'Specific_heat_a={a}_K={K}.jpg')

plt.show()
