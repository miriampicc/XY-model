import numpy as np
from statsmodels.tsa.stattools import acf
import matplotlib.pyplot as plt

# Define the list of temperatures
temperatures = [ 2.0,1.2, 1.0, 0.009]

# Create an empty list to store autocorrelation results
autocorrelations = []

# Load the magnetization data and calculate autocorrelations for each temperature
for temp in temperatures:
    file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output/T_{temp}/Magnetization.txt"
    magnetization_data = np.loadtxt(file_path)
    autocorrelation = acf(magnetization_data, nlags=300)
    autocorrelations.append(autocorrelation)

# Create lags for the x-axis
lags = np.arange(len(autocorrelations[0]))

# Plot all autocorrelation functions on the same plot
plt.figure(figsize=(10, 6))
for i, autocorr in enumerate(autocorrelations):
    label = f'T ={temperatures[i]} K'
    plt.plot(lags, autocorr, linestyle='-', label=label) #marker='o', 

plt.title('Autocorrelation Functions for Different Temperatures')
plt.xlabel('Lag')
plt.ylabel('Autocorrelation')
plt.grid(True)
plt.legend()
plt.show()
