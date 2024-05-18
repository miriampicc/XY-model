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



# Function to calculate the Binder cumulant
def binder_cumulant(m):
    m2 = np.mean(m**2)
    m4 = np.mean(m**4)
    return m4 / (3 * m2**2)

def jackknife_analysis (magns, temp):
    M = np.array(magns)
    n = len(M)
    b_c = np.zeros(n)

    for i in range(n):

        sample = np.delete (M, i)

        b_c [i] = binder_cumulant(sample)

    mean_bc = np.mean(b_c)
    jackknife_error = np.sqrt(((n - 1) / n) * np.sum((b_c - mean_bc)**2))

    return mean_bc, jackknife_error


i=0 

temperatures = []
delta = (1/beta_low - 1/beta_high)/(rank)
T_high = 1/beta_low
T_low = 1/beta_high

mean_values_U = []
jackknife_bias = []
jackknife_errors = [] 
jd_error = []

for l in L : 

    N = l * l

    mean_values_U = []
    jackknife_bias = []
    jackknife_errors = [] 
    jd_error = []
    temperatures = []

    for n in range(rank): 

        magn = []

        t = T_high - n * delta
        print (t)
        
        temperatures.append(t)

        file_path = f"/home/x_mirpi/Output_TBG/K_{K}/e_{e}/L{l}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}/beta_{n}" + '/trsb_magnetization.txt'
        with open(file_path, "r") as file:
            magn = [float(line.strip()) for line in file.readlines()]
    
        n= len (magn)
        jackknife_mean_temp, jackknife_error_temp = jackknife_analysis(magn, t) 
        mean_values_U.append(jackknife_mean_temp)
        jackknife_errors.append(jackknife_error_temp)

        i += 1


    print(f"Processed temperature: {t}")

    mean_val_U = np.array(mean_values_U)
    error_U = np.array(jackknife_errors) 

    plt.plot(temperatures, mean_val_U, label=f'L={l}')
    plt.fill_between(temperatures, mean_val_U - error_U, mean_val_U + error_U, alpha=0.6, linewidth=4)   
            



plt.ylabel(f'Binder Cumulant ')
plt.xlabel('T(K) ')
plt.grid()
plt.legend()
plt.title(f'Binder Cumulant with Jackknife resampling method e={e}, K={K} ')
plt.savefig(f'Binder_Cumulant_JK_e{e}_K{K.jpg}')

plt.show()
