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
temperatures = []
delta = (1/beta_low - 1/beta_high)/(rank)
T_high = 1/beta_low
T_low = 1/beta_high
delta_beta = (beta_high - beta_low)/(rank)

#print(delta)

k=0

for l in L:
    
    N= l * l

    Psi1_r = []
    Psi2_r = []

    mean_Psi1 = [] 
    mean_Psi2 = [] 
    sum_Psis = []

    temperatures = []

    betas = []

    for n in range(rank) : 

        t = T_high - n * delta
        print (t)

        bb = beta_low + delta_beta * n 
        
        temperatures.append(t)
        betas.append(bb)

        Psi1_r = []
        Psi2_r = []


        file_path1 = f"/home/x_mirpi/Output_TBG/K_{K}_tdf/e_{e}/L{l}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}/beta_{n}" + '/Lattice_fin.txt'
        with open(file_path1, "r") as file:
            for line in file:
                columns = line.strip().split()
                if len(columns) == 4:
                    Psi1_r.append(float(columns[0]))
                    Psi2_r.append(float(columns[2]))

        squared_Psi1 = [x**2 for x in Psi1_r]
        squared_Psi2 = [x**2 for x in Psi2_r]

        mm_Psi1 = calculate_mean(squared_Psi1)
        mm_Psi1 = mm_Psi1 / N
        std_psi1 = calculate_std(squared_Psi1)
        mm_Psi2 = calculate_mean(squared_Psi2)
        mm_Psi2 = mm_Psi2 / N
        std_psi2 = calculate_std(squared_Psi2)

        combination = mm_Psi1+ mm_Psi2

        mean_Psi1.append(mm_Psi1)
        mean_Psi2.append(mm_Psi2)
        sum_Psis.append(combination)
        
    Psi1_matrix = np.column_stack((temperatures, mean_Psi1))
    Psi2_matrix = np.column_stack((temperatures, mean_Psi2))
    Sum_matrix = np.column_stack((temperatures, sum_Psis))
    
    plt.plot(Psi1_matrix[:, 0], Psi1_matrix[:, 1], label=r'$|\psi_1|^2$', color = colors[k])
    plt.plot(Psi2_matrix[:, 0], Psi2_matrix[:, 1], label=r'$|\psi_2|^2$', color = colors[k+1])
    plt.plot(Sum_matrix[:, 0], Sum_matrix[:, 1], label=r'$|\psi_2|^2 + |\psi_2|^2$', color = colors[k+2])


    k = k+1

plt.xlabel(r'$T$')
plt.ylabel('Probability density of the compounds')
#plt.title(f'K={K} e = {e} L_0 = {L0}')
plt.legend()
plt.grid(True)
plt.savefig(f'Probability_densities_Psi_tot_L={l}.jpg')

plt.show()