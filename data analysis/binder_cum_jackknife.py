import numpy as np
import matplotlib.pyplot as plt
import statistics

L = [8, 12, 16, 20, 24, 32, 40, 48]
K = 0.5


temperatures = [ 0.5, 0.52, 0.55, 
                0.57, 0.6, 0.62, 0.65, 0.67, 0.7, 0.72, 0.75, 0.77, 0.8, 0.82, 0.85, 0.87, 0.9, 
                0.92, 0.94, 0.95, 0.96, 0.98, 1.0, 1.02, 1.04, 1.05, 1.06, 1.07, 1.09, 1.1, 1.13, 1.15, 1.17, 1.2, 1.22, 1.25, 
                1.27]  

# 1.3, 1.4, 1.45, 1.5, 1.7, 2.0, 2.5
#0.009, 0.01, 0.02, 0.03, 0.05, 0.07, 0.08, 0.09, 0.1, 0.11, 0.13, 0.14, 0.15, 0.17, 0.2, 0.22, 0.25, 0.27, 0.3, 0.32, 0.35, 0.37, 0.4, 0.45, 



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



    for temp in temperatures: 

        magn = []

        file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output_K={K}/Output_L={l}/T_{temp}" + '/trsb_magnetization.txt'
        with open(file_path, "r") as file:
            magn = [float(line.strip()) for line in file.readlines()]
    
        n= len (magn)
        jackknife_mean_temp, jackknife_error_temp = jackknife_analysis(magn, temp) 
        mean_values_U.append(jackknife_mean_temp)
        jackknife_errors.append(jackknife_error_temp)

        i += 1


    print(f"Processed temperature: {temp}")

    mean_val_U = np.array(mean_values_U)
    error_U = np.array(jackknife_errors) 

    plt.plot(temperatures, mean_val_U, label=f'L={l}')
    plt.fill_between(temperatures, mean_val_U - error_U, mean_val_U + error_U, alpha=0.6, linewidth=4)   
            



plt.ylabel(f'Binder Cumulant ')
plt.xlabel('T(K) ')
plt.grid()
plt.legend()
plt.title(f'Binder Cumulant with Jackknife resampling method ')

plt.show()
