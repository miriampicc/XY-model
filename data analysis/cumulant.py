import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 

temperatures = [ 0.009, 0.01, 0.02, 0.03, 0.05, 0.07, 0.08, 0.09, 0.1, 0.11, 0.13, 0.14, 0.15, 
                0.17, 0.2, 0.22, 0.25, 0.27, 0.3, 0.32, 0.35, 0.37, 0.4, 0.45, 0.5, 0.52, 0.55, 
                0.57, 0.6, 0.62, 0.65, 0.67, 0.7, 0.72, 0.75, 0.77, 0.8, 0.82, 0.85, 0.87, 0.9, 
                0.92, 0.94, 0.95, 0.96, 0.98, 1.0, 1.02, 1.04, 1.05, 1.06, 1.07, 1.09, 1.1, 1.13, 
                1.15, 1.17, 1.2, 1.22, 1.25, 1.27, 1.3, 1.4, 1.45, 1.5, 1.7, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5]  

L=[8, 12, 16, 20, 24, 32]
K=1

def calculate_std(data): 
    try:
        std_deviation = np.std(data)
        return std_deviation
    except Exception as e:
        print(f"Error: {e}")
        return None

def calculate_mean(data): 
    mean = sum(data) / len(data)
    return mean


pseudo_magn = []
mm_2_values = []
mm_4_values = [] 
cumulant = []
i = 0
for l in L : 

    N = l * l

    pseudo_magn = []
    mm_2_values = []
    mm_4_values = [] 
    cumulant = []

    for temp in temperatures:
        
        file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Relative_fluctuations/K={K}/Output_L={l}/T_{temp}" + '/trsb_magnetization.txt'

        #file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output_K={K}/Output_L={L}/T_{temp}" + '/trsb_magnetization.txt'
        #file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output_prova/Output_L={L}/T_{temp}" + '/trsb_magnetization.txt'


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

        i += 1
        print(mm)
        print (temp)


    plt.plot(temperatures, pseudo_magn, linestyle='-')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Pseudo- Magnetisation ')
    plt.title(f'Relative density Fluctuations, L={L}')
    plt.grid(True)

plt.show()



pseudo_magn = []
mm_2_values = []
mm_4_values = [] 
cumulant = []
i = 0
for l in L : 

    N = l * l

    pseudo_magn = []
    mm_2_values = []
    mm_4_values = [] 
    cumulant = []

    for temp in temperatures:
        
        file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Relative_fluctuations/K={K}/Output_L={l}/T_{temp}" + '/trsb_magnetization.txt'

        #file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output_K={K}/Output_L={L}/T_{temp}" + '/trsb_magnetization.txt'
        #file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output_prova/Output_L={L}/T_{temp}" + '/trsb_magnetization.txt'


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

        i += 1
        print(mm)
        print (temp)

    plt.plot(temperatures, cumulant, linestyle='-')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Binder cumulant')
    plt.title(f'Relative density Fluctuations, L={L}')
    plt.grid(True)

plt.show()
