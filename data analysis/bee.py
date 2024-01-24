import numpy as np
from statsmodels.tsa.stattools import acf
import matplotlib.pyplot as plt

# Define the list of temperatures
temperatures = [10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.5, 4.0, 3.5, 3.0, 2.5, 2.4, 2.3, 2.2, 2.1, 2.0, 1.9, 1.8, 
                1.7, 1.6, 1.5, 1.45, 1.4, 1.35, 1.3, 1.25, 1.2, 1.15, 1.1, 1.05, 1.0, 0.95, 0.9, 0.85, 0.8, 
                0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.14, 0.13, 0.11, 
                0.1, 0.09, 0.08, 0.07, 0.05, 0.03, 0.02, 0.01, 0.009]
                
               

# Create empty lists to store autocorrelation results and correlation times
autocorrelations = []
correlation_times = []

# Load the magnetization data and calculate autocorrelations for each temperature
for temp in temperatures:
    file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output/T_{temp}/Magnetization.txt"
    magnetization_data = np.loadtxt(file_path)
    autocorrelation = acf(magnetization_data, nlags=50)
    autocorrelations.append(autocorrelation)
    
    # Calculate \tau by finding the lag where autocorrelation drops below a threshold
    threshold = 1 / 2.71828  # 1/e
    tau = None
    for lag, acf_value in enumerate(autocorrelation):
        if acf_value < threshold:
            tau = lag
            break
    correlation_times.append(tau)

lags = np.arange(len(autocorrelations[0]))

# Plot all autocorrelation functions on the same plot
plt.figure(figsize=(10, 6))
for i, autocorr in enumerate(autocorrelations):
    label = f' {temperatures[i]}'
    plt.plot(lags, autocorr, linestyle='-', label=label) #marker='o',

plt.title('Autocorrelation Functions for Different Temperatures')
plt.xlabel('Lag')
plt.ylabel('Autocorrelation')
plt.grid(True)
plt.legend()
plt.show()

# Plot the correlation times as a function of temperature
plt.figure(figsize=(10, 6))
plt.plot(temperatures, correlation_times, linestyle='-') #marker='o', 
plt.title('Correlation Time (\u03C4) as a Function of Temperature')
plt.xlabel('Temperature')
plt.ylabel('Correlation Time (\u03C4)')
plt.grid(True)
plt.show()
