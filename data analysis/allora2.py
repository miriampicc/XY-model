import numpy as np
import matplotlib.pyplot as plt

# List of temperatures
temperatures = [0.01, 0.25, 0.5, 0.55, 0.7, 1.0]

# Load energy values for each temperature (you'll need to adjust the file paths)
energy_values = []

for temp in temperatures:
    energy_values.append(np.loadtxt(f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output/T_{temp}/Energy.txt"))

    steps = len(energy_values)

    # Define the number of bins
    num_bins = 20

    log_bins = np.logspace(0, np.log10(steps), num_bins, base=10.0)

    # Initialize arrays to store mean, std dev, and half width of the bin for each bin
    mean_values = []
    std_dev_values = []
    half_widths = []

    for j in range(len(log_bins) - 1):
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


    plt.scatter(half_widths, mean_values, label=f'T={temp}', marker='o')
    plt.errorbar(half_widths, mean_values, yerr=std_dev_values,linestyle='None', color='black', capsize=3)

plt.xscale('log')  # Set the x-axis to be logarithmic
plt.xlabel('Monte Carlo Steps (log scale)')
plt.ylabel('Values')
plt.legend()
plt.title('Mean Energy Thermalization')
plt.grid()
plt.show()
