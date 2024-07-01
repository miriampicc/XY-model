import numpy as np
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
    mean = np.mean(data)
    return mean

colors = [
    "#00008B",  # Dark Blue
    "#0000CD",  # Medium Blue
    "#1E90FF",  # Dodger Blue
    "#87CEEB",  # Sky Blue
    "#EEE8AA",  # Pale Goldenrod
    "#F0E68C",  # Khaki
    "#FFD700",  # Gold
    "#F08080",  # Light Coral
    "#FF6347"   # Tomato
]

# Let us obtain the temperatures
temperatures = []
delta = (1/beta_low - 1/beta_high)/(rank)
T_high = 1/beta_low
T_low = 1/beta_high

i = 0

for l in L:
    pseudo_magn = []
    cumulant = []
    temperatures = []

    for n in range(rank):
        t = T_high - n * delta
        temperatures.append(t)

        file_path = f"/home/x_mirpi/Output_TBG/K_{K}/e_{e}/L{l}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}/beta_{n}/trsb_magnetization.txt"

        try:
            with open(file_path, 'r') as file:
                numbers = []
                for line_number, line in enumerate(file, start=1):
                    stripped_line = line.strip()
                    try:
                        number = float(stripped_line)
                        numbers.append(number)
                    except ValueError:
                        print(f"Warning: Could not convert line {line_number} to float: {stripped_line}")

                if numbers:
                    m2 = [x**2 for x in numbers]
                    m4 = [x**4 for x in numbers]
                    avg_m2 = calculate_mean(m2)
                    avg_m4 = calculate_mean(m4)
                    mm = calculate_mean(numbers)
                    U = avg_m4 / (3 * avg_m2**2)
                    cumulant.append(np.abs(U))
                    pseudo_magn.append(np.abs(mm))

        except FileNotFoundError:
            print(f"Error: File {file_path} not found.")

    plt.plot(temperatures, cumulant, linestyle='-', label=f'L={l}', color=colors[i])
    i += 1

plt.xlabel('Temperature')
plt.ylabel('Binder Cumulant (U)')
plt.title(f'Binder Cumulant for e={e}')
plt.legend()
plt.grid(True)
plt.savefig(f'Binder_Cumulant_e={e}_K={K}.jpg')
plt.show()
