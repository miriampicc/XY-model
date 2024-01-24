import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 

numbers = ['1', '2']
arguments = ['/Magnetization1.txt', '/Magnetization2.txt'] 
title = ['Lattice 1', 'Lattice 2']

temperatures = [ 0.009, 0.01, 0.02, 0.03, 0.05, 0.07, 0.08, 0.09, 0.1, 0.11, 0.13, 0.14, 0.15, 
                0.17, 0.2, 0.22, 0.25, 0.27, 0.3, 0.32, 0.35, 0.37, 0.4, 0.45, 0.5, 0.52, 0.55, 
                0.57, 0.6, 0.62, 0.65, 0.67, 0.7, 0.72, 0.75, 0.77, 0.8, 0.82, 0.85, 0.87, 0.9, 
                0.92, 0.94, 0.95, 0.96, 0.98, 1.0, 1.02, 1.04, 1.05, 1.06, 1.07, 1.09, 1.1, 1.13, 1.15, 1.17, 1.2, 1.22, 1.25, 
                1.27, 1.3, 1.4, 1.45, 1.5, 1.7, 2.0, 2.5]  

N = 12 * 12

# Create a dictionary to store the arrays
mean_argument_values = []

def calculate_std(data): 
    try:
        std_deviation = np.std(data)
        return std_deviation
    except Exception as e:
        # Handle any exceptions that may occur during the calculation
        print(f"Error: {e}")
        return None

def calculate_mean(data): 
    mean = sum(data) / len(data)
    return mean



#PLOT OF THE ENERGY 
i = 0
for temp in temperatures:
        

    file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output/T_{temp}" + '/Energy.txt'

    with open(file_path, 'r') as file:
        numbers = [float(line.strip()) for line in file.readlines()]
        mm = calculate_mean(numbers)
        mean_argument_values.append(mm) 

    i += 1
    #print(i)


np.savetxt ('/Users/mirimi/Desktop/hihi/KTH/XY-model/data analysis/n=100000/J1=2_J2=K/L=12/Energy.txt', mean_argument_values,  header='Energy', comments='')

# Plot Energy vs. Temperature
plt.plot(temperatures, mean_argument_values, marker='o', linestyle='-')
plt.xlabel('Temperature (K)')
plt.ylabel('Energy')
plt.title('Mean Energy  ($J_1 =2 J_2=K$), n$_{steps}$= 100000')
plt.grid(True)

plt.show()


# Create a dictionary to store the arrays
mean_argument_values = {argument: [] for argument in arguments}

#PLOT OF THE MAGNETISATION FIRST LAYER
i = 0
for argument in arguments: 

    mean_argument_values[argument] = []

    print(argument)

    for temp in temperatures:
        

        file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output/T_{temp}" + argument

        with open(file_path, 'r') as file:
            numbers = [float(line.strip()) for line in file.readlines()]
            mm = calculate_mean(numbers)
            mean_argument_values[argument].append(mm) 

        i += 1
        #print(i)
    
    np.savetxt ('/Users/mirimi/Desktop/hihi/KTH/XY-model/data analysis/n=100000/J1=2_J2=K/L=12' + argument , mean_argument_values[argument],  header='Magnetisation', comments='')

   
i=0
for argument, values in mean_argument_values.items():
    plt.plot(temperatures, values, marker='o', linestyle='-', label= title[i] )
    i += 1

plt.xlabel('Temperature (K)')
plt.ylabel('Magnetisation')
plt.title('Mean Magnetisation ($J_1 =2 J_2=K$), n$_{steps}$= 100000')
plt.grid(True)
plt.legend()
plt.show()


#PLOT OF TRSB MAGNETISATION
# Create a dictionary to store the arrays
mean_argument_values = []
mm_2_values = []
mm_4_values = [] 
cumulant = []
i = 0
for temp in temperatures:
        

    file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output/T_{temp}" + '/trsb_magnetization.txt'

    with open(file_path, 'r') as file:
        numbers = [float(line.strip()) for line in file.readlines()]
        m2= [x**2 for x in numbers ]
        m4= [x**4 for x in numbers ]
        avg_m2 = calculate_mean (m2)
        avg_m4 = calculate_mean (m4)
        mm = calculate_mean(numbers)
        U = avg_m4 / (3 * avg_m2**2)
        cumulant.append (np.abs(U))
        mean_argument_values.append(np.abs(mm))

    i += 1
    print(mm)
    print (temp)

np.savetxt('/Users/mirimi/Desktop/hihi/KTH/XY-model/data analysis/n=100000/J1=2_J2=K/L=12/cumulant.txt', cumulant, header='Cumulant values', comments='')
np.savetxt('/Users/mirimi/Desktop/hihi/KTH/XY-model/data analysis/n=100000/J1=2_J2=K/L=12/pseudo_magn.txt', mean_argument_values, header='Pseudo Magnetisation values', comments='')


# Plot Energy vs. Temperature
plt.plot(temperatures, mean_argument_values, marker='o', linestyle='-')
plt.xlabel('Temperature (K)')
plt.ylabel('Pseudo- Magnetisation ')
plt.title('Pseudo- Magnetisation   ($J_1 =2 J_2=K$), n$_{steps}$= 100000')
plt.grid(True)

plt.show()

#Cumulant Binder 
plt.plot(temperatures, cumulant, marker='o', linestyle='-')
plt.xlabel('Temperature (K)')
plt.ylabel('Binder cumulant')
plt.title('Binder Cumulant   ($J_1 =2 J_2=K$), n$_{steps}$= 100000')
plt.grid(True)

plt.show()


# Superfluid Stiffness
Jp1 = []
Js1 = []
mean_Jd_values1 = [] 
sin_Ic1= []

i = 0 

for temp in temperatures:
    Jd1 = []
    Ic1 = []

    file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output/T_{temp}" + '/Helicity_modulus1.txt'
    with open(file_path, "r") as file:
        for line in file:
            columns = line.strip().split()
            if len(columns) == 2:
                Jd1.append(float(columns[0]))
                Ic1.append(float(columns[1]))

    mean_sin = calculate_std(Ic1) ** 2
    sin_Ic1.append(mean_sin)
    sin = mean_sin * N / temperatures[i]

    Jp1.append(sin)

    mm = calculate_mean(Jd1) 
    cos_Jd = mm 
    mean_Jd_values1.append(cos_Jd)

    Js_new = cos_Jd - sin
    Js1.append(Js_new)
    
    i += 1 

np.savetxt('/Users/mirimi/Desktop/hihi/KTH/XY-model/data analysis/n=100000/J1=2_J2=K/L=12/Jp1.txt', Jp1, header='Paramagnetic contribution1', comments='')
np.savetxt('/Users/mirimi/Desktop/hihi/KTH/XY-model/data analysis/n=100000/J1=2_J2=K/L=12/Jd1.txt', mean_Jd_values1, header='Diamagnetic contribution1', comments='')
np.savetxt('/Users/mirimi/Desktop/hihi/KTH/XY-model/data analysis/n=100000/J1=2_J2=K/L=12/Js1.txt', Js1, header='Superfluidity stiffness1', comments='')


Jd_matrix1 = np.column_stack((temperatures, mean_Jd_values1))
Jp_matrix1 = np.column_stack((temperatures, Jp1))
Js_matrix1 = np.column_stack((temperatures, Js1))

y = [0] * len(temperatures)

for temp in range(len(temperatures)): 
    y[temp] = 2 * temperatures[temp] / np.pi

plt.plot(Jd_matrix1[:, 0], Jd_matrix1[:, 1], linestyle='-', label=r'$J_d$')
plt.plot(Jp_matrix1[:, 0], Jp_matrix1[:, 1],  linestyle='-', label = r'$J_p$')
plt.plot(Js_matrix1[:, 0], Js_matrix1[:, 1], linestyle='-', label = r'$J_s$')
plt.plot(temperatures , y, linestyle='--', label = r'$\frac{T}{\pi}$')
plt.legend()
plt.xlabel('Temperature (K)')
plt.ylabel('superfluid stiffness')
plt.title('Superfulitity Stiffness Lattice 1, ($J_1 =2 J_2=K$), n$_{steps}$= 100000')
plt.grid(True)

plt.show()


# Superfluid Stiffness
Jp2 = []
Js2 = []
mean_Jd_values2 = [] 
sin_Ic2= []

i = 0 

for temp in temperatures:
    Jd2 = []
    Ic2 = []

    file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output/T_{temp}" + '/Helicity_modulus2.txt'
    with open(file_path, "r") as file:
        for line in file:
            columns = line.strip().split()
            if len(columns) == 2:
                Jd2.append(float(columns[0]))
                Ic2.append(float(columns[1]))

    mean_sin = calculate_std(Ic2) ** 2
    sin_Ic2.append(mean_sin)
    sin = mean_sin * N / temperatures[i]

    Jp2.append(sin)

    mm = calculate_mean(Jd2) 
    cos_Jd = mm 
    mean_Jd_values2.append(cos_Jd)

    Js_new = cos_Jd - sin
    Js2.append(Js_new)
    
    i += 1 

np.savetxt('/Users/mirimi/Desktop/hihi/KTH/XY-model/data analysis/n=100000/J1=2_J2=K/L=12/Jp2.txt', Jp2, header='Paramagnetic contribution2', comments='')
np.savetxt('/Users/mirimi/Desktop/hihi/KTH/XY-model/data analysis/n=100000/J1=2_J2=K/L=12/Jd2.txt', mean_Jd_values2, header='Diamagnetic contribution2', comments='')
np.savetxt('/Users/mirimi/Desktop/hihi/KTH/XY-model/data analysis/n=100000/J1=2_J2=K/L=12/Js2.txt', Js2, header='Superfluidity stiffness2', comments='')

Jd_matrix2 = np.column_stack((temperatures, mean_Jd_values2))
Jp_matrix2 = np.column_stack((temperatures, Jp2))
Js_matrix2 = np.column_stack((temperatures, Js2))

y = [0] * len(temperatures)

for temp in range(len(temperatures)): 
    y[temp] = 2* temperatures[temp] / np.pi

plt.plot(Jd_matrix2[:, 0], Jd_matrix2[:, 1], linestyle='-', label=r'$J_d$')
plt.plot(Jp_matrix2[:, 0], Jp_matrix2[:, 1],  linestyle='-', label = r'$J_p$')
plt.plot(Js_matrix2[:, 0], Js_matrix2[:, 1], linestyle='-', label = r'$J_s$')
plt.plot(temperatures , y, linestyle='--', label = r'$\frac{T}{\pi}$')
plt.legend()
plt.xlabel('Temperature (K)')
plt.ylabel('superfluid stiffness')
plt.title('Superfulitity Stiffness Lattice 2, ($J_1 =2 J_2=K$), n$_{steps}$= 100000 ')
plt.grid(True)

plt.show()








