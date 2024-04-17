import numpy as np
import matplotlib.pyplot as plt

temperatures =  [ 0.01, 0.25, 0.5 , 0.55 , 0.7, 1.0 ]
#0.009 0.01, 0.02, 0.03, 0.05, 0.07, 0.08, 0.09, 0.1, 0.11, 0.13, 0.14, 0.15, 0.2, 0.25, 0.3, 0.35, 
 #               0.4 , 0.45 , 0.5 , 0.55 , 0.6 , 0.65 , 0.7 , 0.75 , 0.8 , 0.85 , 0.9 , 0.92 , 0.94 , 0.95 , 0.96 , 0.98 , 1.0 , 1.02 , 1.04 , 1.05 , 1.07 ,
  #               1.09 , 1.1 , 1.12 , 1.14 , 1.15 , 1.17 , 1.18 , 1.2 , 1.23 , 1.25 , 1.27 , 1.29 , 
   #             1.3 , 1.33 , 1.35 , 1.37 , 1.4 , 1.45 , 1.5 , 1.6 , 1.7 , 1.8 , 1.9 , 2.0 , 2.1 , 2.2 , 2.3 , 2.4 , 2.5 , 
    #            3.0 , 3.5 , 4.0 , 4.5 , 5.0

energy_values = np.loadtxt("/Users/mirimi/Desktop/hihi/KTH/XY-model/Output/T_0.01/Energy.txt") 
print(len(energy_values))
steps = len(energy_values)

# Define the number of bins
num_bins = 20

log_bins = np.logspace(0, np.log10(steps), num_bins, base=10.0) 

# Initialize arrays to store mean, std dev, and half width of the bin for each bin
mean_values = []
std_dev_values = []
half_widths = []
inside_bin = []
x_mc_val = []

i = 0 
j = 0


for j in range(len(log_bins)-1) : 
    inside_bin = []

    for i in range(len(energy_values)):
        if log_bins[j] <= i < log_bins[j + 1]:
            inside_bin.append(energy_values[i])

    mean_val = np.mean(inside_bin)
    std_val = np.std(inside_bin)
    x_middle = (log_bins[j + 1] + log_bins[j]) / 2.0

    mean_values.append(mean_val)
    std_dev_values.append(std_val)
    half_widths.append(x_middle)

print (half_widths)

plt.scatter(half_widths, mean_values, label='Mean Energy ', marker='o')
plt.errorbar(half_widths, mean_values, yerr=std_dev_values,linestyle='None', color='black', capsize=3, label='Error bars')
plt.xscale('log')  # Set the x-axis to be logarithmic
plt.xlabel('Monte Carlo Steps (log scale)')
plt.ylabel('Energy')
plt.legend()
plt.title('Energy Thermalization, T=0.35, n$_{steps} =1000000$ ')
plt.grid()
plt.show()
     
