import numpy as np
import matplotlib.pyplot as plt
import statistics

L=16

N = L * L
def calculate_mean(data):
    return np.mean(data)

def calculate_standard_deviation(data):
    try:
        stdev = statistics.stdev(data)
        return stdev
    except statistics.StatisticsError:
        return "Input data should contain at least two data points for standard deviation calculation."
    

def jackknife_analysis(Ics, temp):

    Ic = np.array (Ics)

    n = len(Ic)
    Ic_std = (calculate_standard_deviation(Ic))**2

    jackknife_Jp_values = []

    for i in range(n):
        
        data_subset = np.delete(Ic, i)
        #data_subset_mean = np.mean(data_subset)


        squared_std_subset = ((np.std(data_subset))**2) * (N /temp)
        jackknife_Jp_values.append(squared_std_subset)

    jackknife_mean = np.mean(jackknife_Jp_values)
    jackknife_error = np.sqrt(((n - 1) / n) * np.sum((jackknife_Jp_values - jackknife_mean)**2))

    print(len(data_subset))
    #print(len(jackknife_error))

    return jackknife_mean, jackknife_error

    

temperatures = [ 0.009, 0.01, 0.02, 0.03, 0.05, 0.07, 0.08, 0.09, 0.1, 0.11, 0.13, 0.14, 0.15, 
                0.17, 0.2, 0.22, 0.25, 0.27, 0.3, 0.32, 0.35, 0.37, 0.4, 0.45, 0.5, 0.52, 0.55, 
                0.57, 0.6, 0.62, 0.65, 0.67, 0.7, 0.72, 0.75, 0.77, 0.8, 0.82, 0.85, 0.87, 0.9, 
                0.92, 0.94, 0.95, 0.96, 0.98, 1.0, 1.02, 1.04, 1.05, 1.06, 1.07, 1.09, 1.1, 1.13, 1.15, 1.17, 1.2, 1.22, 1.25, 
                1.27, 1.3, 1.4, 1.45, 1.5, 1.7, 2.0, 2.5]  


mean_Jd_values = []

i=0 



mean_values_Jp = []
mean_values_Jd = []
jackknife_bias = []
jackknife_errors = [] 
jd_error = []

for temp in temperatures:

    Jd = []
    Ic = []

    file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output_L={L}/T_{temp}" + '/Helicity_modulus1.txt' 
    with open(file_path, "r") as file:
        for line in file:
            
            columns = line.strip().split()
            if len(columns) == 2:
                
                Jd.append(float(columns[0]))
                Ic.append(float(columns[1]))
    
    n = len(Ic)

    jackknife_mean_temp, jackknife_error_temp = jackknife_analysis(Ic, temp)
    mean_values_Jp.append(jackknife_mean_temp)
    jackknife_errors.append(jackknife_error_temp)

    mean_Jd = calculate_mean(Jd) #
    std_Jd = calculate_standard_deviation (Jd) / (np.sqrt(N-1))

    mean_values_Jd.append(mean_Jd)
    jd_error.append(std_Jd)

    i += 1


    print(f"Processed temperature: {temp}")



mean_val_Jp = np.array(mean_values_Jp)
mean_val_Jd = np.array(mean_values_Jd)
Js = mean_val_Jd - mean_val_Jp

errors_array = np.array (jackknife_errors)
jd_error_array = np.array (jd_error)
error_Js = np.sqrt( jd_error_array **2 + errors_array**2 )


plt.figure()
plt.plot(temperatures, mean_val_Jp, label='Jp')
plt.plot(temperatures, mean_val_Jd, label='Jd')
plt.plot(temperatures, Js, label='Js')
plt.fill_between(temperatures, Js - error_Js, Js + error_Js, alpha=0.6, facecolor='lightgreen', linewidth=4, label='Error Bounds ')
plt.fill_between(temperatures,  mean_val_Jp - errors_array,  mean_val_Jp + errors_array, alpha=0.6, facecolor='lightblue', linewidth=4, label='Error Bounds ')
plt.ylabel(f'Supefulidity Stiffness components ')
plt.xlabel('T(K) ')
plt.grid()
plt.legend()
plt.title(f'Supefulidity Stiffness Lattice 2 with Jackknife resampling method ')

plt.show()
