import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 

temperatures = [ 0.009, 0.01, 0.02, 0.03, 0.05, 0.07, 0.08, 0.09, 0.1, 0.11, 0.13, 0.14, 0.15, 
                0.17, 0.2, 0.22, 0.25, 0.27, 0.3, 0.32, 0.35, 0.37, 0.4, 0.45, 0.5, 0.52, 0.55, 
                0.57, 0.6, 0.62, 0.65, 0.67, 0.7, 0.72, 0.75, 0.77, 0.8, 0.82, 0.85, 0.87, 0.9, 
                0.92, 0.94, 0.95, 0.96, 0.98, 1.0, 1.02, 1.04, 1.05, 1.06, 1.07, 1.09, 1.1, 1.13, 1.15, 1.17, 1.2, 1.22, 1.25, 
                1.27, 1.3, 1.4, 1.45, 1.5, 1.7, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5]  

#L = [8, 12]
L=[8, 12, 16, 20]
#L=40


#K=[11]
K=1

def calculate_mean(data): 
    mean = sum(data) / len(data)
    return mean



#PLOT OF THE SPECIFIC HEAT 
i = 0
for l in L:
    N = l * l

    specific_heat = []
    beta = []

    for temp in temperatures:
        

        #file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output_K=11/Output_L={l}/T_{temp}" + '/Energy.txt'
        file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Relative_fluctuations/K={K}/Output_L={l}/T_{temp}" + '/Energy.txt'


        #file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output_prova/Output_L={l}/T_{temp}" + '/Energy.txt'

        with open(file_path, 'r') as file:
            numbers = [float(line.strip()) for line in file.readlines()]
            mm = calculate_mean(numbers)
            en_var = (np.std(numbers))**2
            cc = ((1/temp)**2) * en_var * (1/N)
        
            specific_heat.append(cc)
        
        
        beta.append (temp)
    
    sh_val = np.array (specific_heat)
    beta_val = np.array (beta)

    # Plot Energy vs. Temperature
    plt.plot(beta_val, sh_val, linestyle='-', label = f'L={l}')  

""" 
k=11

N = L * L

specific_heat = []
beta = []

for temp in temperatures:
        

    file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output_K={k}/Output_L={L}/T_{temp}" + '/Energy.txt'

    with open(file_path, 'r') as file:
        numbers = [float(line.strip()) for line in file.readlines()]
        mm = calculate_mean(numbers)
        en_var = (np.std(numbers))**2
        cc = ((1/temp)**2) * en_var
        
        specific_heat.append(cc)
        
        
    beta.append (temp)
    
sh_val = np.array (specific_heat)
beta_val = np.array (beta)

    # Plot Energy vs. Temperature
plt.plot(beta_val, sh_val, linestyle='-', label = f'$J_1=0.5$')  

 """
   

plt.xlabel(r'$T(K)$')
plt.ylabel('$C_V$')
plt.title(f'$K = {1}$ ')
plt.legend()
plt.grid(True)

plt.show()