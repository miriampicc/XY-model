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
                1.27, 1.3, 1.4, 1.45, 1.5]   
L=[8, 12, 16, 20, 24, 32, 40]
K=1

rainbow_colors = ['red', 'orange', 'yellow', 'green', 'blue', 'indigo', 'violet']


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


def multiple_jackknife (Ics1, Ics2): 

    Ic1 = np.array (Ics1)
    Ic2 = np.array (Ics2)

    mult_12 = []
    molt_mean_subset = []
    mean_molt_subset = []

    for i in range(len(Ic1)): 
        Ic1_subset = np.delete (Ic1, i)
        Ic2_subset = np.delete (Ic2, i)

        for j in range (len(Ic1_subset)): 
            mult_12.append ( Ic1_subset[j] * Ic2_subset[j])

        mean_molt_subset.append (calculate_mean (mult_12))

        molt_mean_subset.append (calculate_mean(Ic1_subset) * calculate_mean (Ic2_subset))
    
    mean_molt_array = np.array (mean_molt_subset)
    molt_mean_array = np.array (molt_mean_subset)
    sub = mean_molt_array - molt_mean_array

    jackknife_mult = calculate_mean (sub)
    jackknife_error = np.sqrt(((n - 1) / n) * np.sum((sub - jackknife_mult)**2))

    return jackknife_mult, jackknife_error


def jackknife_analysis(Ics, temp):

    Ic = np.array (Ics)

    n = len(Ic)
    Ic_std = (calculate_std(Ic))**2

    jackknife_Jp_values = []

    for i in range(n):
        
        data_subset = np.delete(Ic, i)
        #data_subset_mean = np.mean(data_subset)


        squared_std_subset = ((np.std(data_subset))**2) * (N /temp)
        jackknife_Jp_values.append(squared_std_subset)

    jackknife_mean = np.mean(jackknife_Jp_values)
    jackknife_error = np.sqrt(((n - 1) / n) * np.sum((jackknife_Jp_values - jackknife_mean)**2))

    #print(len(data_subset))
    #print(len(jackknife_error))

    return jackknife_mean, jackknife_error


j = 0

for l in L:

    N= l*l

    Jp1 = []
    Js1 = []
    Jd1_values = [] 

    Jp2 = []
    Js2 = []
    Jd2_values = [] 

    double = []
    errors = []

    print(l)

    for temp in temperatures:
        Jd1 = []
        Ic1 = []

        Jd2 = []
        Ic2 = []

    
        file_path1 = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Relative_fluctuations/K={K}/Output_L={l}/T_{temp}" + '/Helicity_modulus1.txt'

        with open(file_path1, "r") as file:
            for line in file:
                columns = line.strip().split()
                if len(columns) == 2:
                    Jd1.append(float(columns[0]))
                    Ic1.append(float(columns[1]))
        

        #Jackknife Ic (error Ic)
        n = len (Ic1)
        jackknife_Jp1, jackknife_error_Jp1 = jackknife_analysis(Ic1, temp)
        print(temp)
        Jp1.append(jackknife_Jp1)

        cos_Jd1 = calculate_mean(Jd1) 
        Jd1_values.append(cos_Jd1)
        error_Jd1 = calculate_std (Jd1) / (np.sqrt(N-1))
        
        Js_new1 = cos_Jd1 - jackknife_Jp1
        Js1.append(Js_new1)
        error_Js1 = np.sqrt(jackknife_error_Jp1**2 + error_Jd1**2)
        

        file_path2 = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Relative_fluctuations/K={K}/Output_L={l}/T_{temp}" + '/Helicity_modulus2.txt'
        #file_path2 = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output_prova/Output_L={L}/T_{temp}" + '/Helicity_modulus2.txt'
        #file_path2 = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output_K={K}/Output_L={L}/T_{temp}" + '/Helicity_modulus2.txt'
        with open(file_path2, "r") as file:
            for line in file:
                columns = line.strip().split()
                if len(columns) == 2:
                    Jd2.append(float(columns[0]))
                    Ic2.append(float(columns[1]))


        n = len (Ic2)
        jackknife_Jp2, jackknife_error_Jp2 = jackknife_analysis(Ic2, temp)
        print(temp)
        Jp2.append(jackknife_Jp2)

        cos_Jd2 = calculate_mean (Jd2)
        Jd2_values.append(cos_Jd2)
        error_Jd2 = calculate_std(Jd2)/ (np.sqrt(N-1)) 

        Js_new2 = cos_Jd2 - jackknife_Jp2
        Js2.append(Js_new2)
        error_Js2 = np.sqrt(jackknife_error_Jp2**2 + error_Jd2**2)


        #Double component
        Js_12 , Js_12_error = multiple_jackknife (Ic1, Ic2)
        print(temp)
        sub = Js_12 
        sott = Js_new1 + Js_new2 +2 * sub * 1/temp
        tot_error = np.sqrt( error_Js1**2 +  error_Js2**2 + Js_12_error**2)


    errors.append (tot_error)
    double.append(sott)
    


    Jd_matrix1 = np.column_stack((temperatures, Jd1_values))
    Jp_matrix1 = np.column_stack((temperatures, Jp1))
    Js_matrix1 = np.column_stack((temperatures, Js1))

    Jd_matrix2 = np.column_stack((temperatures, Jd2_values))
    Jp_matrix2 = np.column_stack((temperatures, Jp2))
    Js_matrix2 = np.column_stack((temperatures, Js2))

    duouble_matrix = np.column_stack ((temperatures, double))
    double_matrix_error = np.column_stack ((temperatures, errors))

    y = [0] * len(temperatures)

    for temp in range(len(temperatures)): 
        y[temp] = 2 * temperatures[temp] / np.pi


    """plt.figure(figsize=(15, 5))

    # Plot 1
    plt.subplot(1, 3, 1)  # 3 rows, 1 column, plot 1
    plt.plot(Jd_matrix1[:, 0], Jd_matrix1[:, 1], label='Jd')
    plt.plot(Jp_matrix1[:, 0], Jp_matrix1[:, 1], label='Jp')
    plt.plot(Js_matrix1[:, 0], Js_matrix1[:, 1], label='Js')
    plt.plot(temperatures , y, linestyle='--', label = r'$\frac{2T}{\pi}$')

    plt.xlabel('Temperature (K)')
    plt.ylabel('superfluid stiffness')
    plt.title(' Lattice 1')
    plt.legend()
    plt.grid(True)

    # Plot 2
    plt.subplot(1, 3, 2)  # 3 rows, 1 column, plot 2
    plt.plot(Jd_matrix2[:, 0], Jd_matrix2[:, 1], label='Jd')
    plt.plot(Jp_matrix2[:, 0], Jp_matrix2[:, 1], label='Jp')
    plt.plot(Js_matrix2[:, 0], Js_matrix2[:, 1], label='Js')
    plt.plot(temperatures , y, linestyle='--', label = r'$\frac{2T}{\pi}$')

    plt.legend()
    plt.xlabel('Temperature (K)')
    plt.ylabel('superfluid stiffness')
    plt.title('Lattice 2')
    plt.grid(True)

    # Plot 3
    plt.subplot(1, 3, 3)  # 3 rows, 1 column, plot 3"""
    plt.plot(duouble_matrix[:, 0], duouble_matrix[:, 1], label = f'L ={l}', color=rainbow_colors[j % len(rainbow_colors)])
    #plt.fill_between (duouble_matrix[:, 0], duouble_matrix[:, 1]- double error )

    j += 1

plt.plot(temperatures , y, linestyle='--', label = r'$\frac{2T}{\pi}$')
plt.legend()
plt.xlabel('Temperature (K)')
plt.ylabel('Helicity sum')
plt.title('Helicity sum of bilayer compound')
plt.grid(True)

# Add a title to the entire figure
#plt.suptitle('Relative density fluctuations K,=1 L=16')

#plt.tight_layout()  # Adjust layout for better appearance



plt.show()
