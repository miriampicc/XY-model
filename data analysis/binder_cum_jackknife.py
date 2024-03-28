import numpy as np
import matplotlib.pyplot as plt
import statistics

L = [8, 12, 16, 20, 24, 32, 40, 64] #, 40, 64
K = 1


temperatures = [  0.2, 0.22, 0.25, 0.27, 0.3, 0.32, 0.35, 0.37, 0.4, 0.45, 0.5, 0.52, 0.55, 
                0.57, 0.6, 0.62, 0.65, 0.67, 0.7, 0.72, 0.75, 0.77, 0.8, 0.82, 0.85, 0.87, 0.9] 

"""temperatures = [0.4, 0.405, 0.41, 0.415, 0.42, 0.422, 0.424, 0.426, 0.428, 0.43, 0.4305, 0.4310, 0.4315, 0.4320, 0.4325,
                0.4330, 0.4335, 0.4340, 0.4345, 0.4350, 0.4355, 0.4360, 0.4365, 0.4370, 0.4375, 0.4380, 0.4385, 0.4390, 0.4395, 0.44,
                0.4405, 0.4410, 0.4415, 0.4420, 0.4425, 0.4430, 0.4435, 0.4440, 0.4445, 0.4450, 0.4455, 0.4460, 0.4465, 0.4470, 0.4475,
                0.4480, 0.4485, 0.4490, 0.4495, 0.45, 0.4505, 0.4510, 0.4515, 0.4520, 0.4525, 0.4530, 0.4535, 0.4540, 0.4545, 0.4550,
                0.4555, 0.4560, 0.4565, 0.4570, 0.4575, 0.4580, 0.4585, 0.4590, 0.4595, 0.46]"""

rainbow_colors = ['red', 'orange', 'yellow', 'green', 'blue', 'indigo', 'violet', 'brown']

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

b=0

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


        file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Relative_fluctuations/K={K}/Output_L={l}/T_{temp}" + '/trsb_magnetization.txt'

        #file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output_K={K}/Output_L={l}/T_{temp}" + '/trsb_magnetization.txt'
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

    plt.plot(temperatures, mean_val_U, label=f'L={l}', color=rainbow_colors[b % len(rainbow_colors)])
    plt.fill_between(temperatures, mean_val_U - error_U, mean_val_U + error_U, alpha=0.6, linewidth=4)  
    b += 1 
            



plt.ylabel(f'Binder Cumulant ')
plt.xlabel('T(K) ')
plt.grid()
plt.legend()
plt.title(f'Binder Cumulant with Jackknife resampling method ')

plt.show()
