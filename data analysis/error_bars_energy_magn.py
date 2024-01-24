import numpy as np
import statistics
import matplotlib.pyplot as plt


arguments = ['/Energy.txt', '/Magnetization1.txt', '/Magnetization2.txt'] 
titles = [ 'Energy', 'Magnetization1', 'Magnetization2', 'trsb magnetization']


def calculate_mean(data):
    return np.mean(data)

def load_energy_data(temp , argument):
    file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output_L=20/T_{temp}" + argument
    with open(file_path, 'r') as file:
        return [float(line.strip()) for line in file.readlines()]
    
    

def calculate_standard_deviation(data):
    try:
        stdev = statistics.stdev(data)
        return stdev
    except statistics.StatisticsError:
        return "Input data should contain at least two data points for standard deviation calculation."



temperatures = [ 0.009, 0.01, 0.02, 0.03, 0.05, 0.07, 0.08, 0.09, 0.1, 0.11, 0.13, 0.14, 0.15, 
                0.17, 0.2, 0.22, 0.25, 0.27, 0.3, 0.32, 0.35, 0.37, 0.4, 0.45, 0.5, 0.52, 0.55, 
                0.57, 0.6, 0.62, 0.65, 0.67, 0.7, 0.72, 0.75, 0.77, 0.8, 0.82, 0.85, 0.87, 0.9, 
                0.92, 0.94, 0.95, 0.96, 0.98, 1.0, 1.02, 1.04, 1.05, 1.06, 1.07, 1.09, 1.1, 1.13, 1.15, 1.17, 1.2, 1.22, 1.25, 
                1.27, 1.3, 1.4, 1.45, 1.5, 1.7, 2.0, 2.5]  

mean_argument_values = {argument: [] for argument in arguments}

i=0
for argument in arguments: 
    mean_argument_values[argument] = []
    mean_values = []
    jackknife_bias = []
    errors = []

    for temp in temperatures:

        energies = load_energy_data(temp, argument)

        N = len(energies)
    
        if energies is not None:

            mm = calculate_mean (energies)
            std = calculate_standard_deviation(energies) / (np.sqrt(N-1))

            mean_values.append(mm)
            errors.append(std)

    mean_val = np.array(mean_values)
    errors_array = np.array (errors)

    plt.figure()
    plt.plot(temperatures, mean_val, label='Mean' f'{titles[i]}')
    plt.fill_between(temperatures, mean_val-errors_array, mean_val+errors_array, alpha=0.2, facecolor='b', linewidth=4, label='Error Bounds')
    plt.xlabel('Temperature (K)')
    plt.ylabel(f'{titles[i]}')
    plt.legend()
    plt.title(f'Mean {titles[i]} vs. ($J_1 = J_2=K=1 $), n= 100000')
    plt.grid(True)
    i += 1

plt.show()

mean_argument_values[argument] = []
mean_values = []
jackknife_bias = []
errors = []

for temp in temperatures:

    trsb = load_energy_data(temp, '/trsb_magnetization.txt')

    N = len(trsb)
    
    if trsb is not None:

        mm = np.abs(calculate_mean (trsb))
        std = np.abs(calculate_standard_deviation(trsb) / (np.sqrt(N-1)))

        mean_values.append(mm)
        errors.append(std)

mean_val = np.array(mean_values)
errors_array = np.array (errors)

plt.figure()
plt.plot(temperatures, mean_val, label='Mean trsb magnetisation')
plt.fill_between(temperatures, mean_val-errors_array, mean_val+errors_array, alpha=0.2, facecolor='b', linewidth=4, label='Error Bounds')
plt.xlabel('Temperature (K)')
plt.ylabel('TRSB magnetisation')
plt.legend()
plt.title('TRSB magnetisation vs. ($J_1 = J_2=K=1 $), n= 100000')
plt.grid(True)

plt.show()
