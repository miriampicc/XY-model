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
#print(delta)

delta_beta = (beta_high - beta_low)/(rank)

i=0

for l in L:
    
    N= l * l

    pseudo_magn = []
    mm_2_values = []
    mm_4_values = [] 
    cumulant = []
    temperatures = []

    betas = []

    for n in range(rank) : 

        t = T_high - n * delta
        print (t)

        bb = beta_low + delta_beta * n 
        
        temperatures.append(t)
        betas.append(bb)

        file_path = f"/home/x_mirpi/Output_TBG/K_{K}/e_{e}/L{l}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}/beta_{n}" + '/trsb_magnetization.txt'

        #print(n)

        with open(file_path, 'r') as file:
            numbers = [float(line.strip()) for line in file.readlines()]
            m2= [x**2 for x in numbers ]
            m4= [x**4 for x in numbers ]
            avg_m2 = calculate_mean (m2)
            avg_m4 = calculate_mean (m4)
            mm = calculate_mean(numbers)
            U = avg_m4 / (3 * avg_m2**2)
            cumulant.append (np.abs(U))
            pseudo_magn.append(np.abs(mm))

    beta_array = np.array(betas)

    
    plt.plot(beta_array, cumulant, linestyle='-', label = f'L={l}', color=colors[i])

    i = i+1


plt.xlabel(r'$\beta$')
plt.ylabel('U')
plt.title(f'Binder Cumulant e={e}')
plt.legend()
plt.grid(True)
plt.savefig(f'Binder_Cumulant_vs_beta_e={e}_K={K}.jpg')


plt.show()