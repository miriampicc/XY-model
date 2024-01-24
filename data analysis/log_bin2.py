import numpy as np
from statsmodels.tsa.stattools import acf
import matplotlib.pyplot as plt

# Define the list of temperatures
temperatures = [2.0, 1.2, 1.0, 0.009]
taus = []

# Create an empty list to store autocorrelation results
autocorrelations = []
threshold = 1 / 2.71828

# Load the magnetization data and calculate autocorrelations for each temperature
for temp in temperatures:
    file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output/T_{temp}/Magnetization.txt"
    magnetization_data = np.loadtxt(file_path)
    autocorrelation = acf(magnetization_data, nlags=250)
    autocorrelations.append(autocorrelation)
    nlags=250 
    
    tau = None
    for lag, acf_value in enumerate(autocorrelation):
        if acf_value < threshold:
            tau = lag
            break
    taus.append(tau)

    print(f'T={temp}K with correlation time (tau)={tau}')  # Print correlation time for each temperature

    #print(f'T={temp}K with $\tau$={tau}')

# Total number of Monte Carlo steps
total_steps = len(magnetization_data)

# Create time in terms of Monte Carlo steps for the x-axis
time = np.arange(len(autocorrelations[0])) * (total_steps / nlags)
tau_time = [tau * (total_steps / nlags) for tau in taus]

for temp, tau_phys in zip(temperatures, tau_time):
    print(f'T={temp}K with correlation time (tau)={tau_phys} MC steps')

ttau = np.array(tau_time)



# Plot all autocorrelation functions on the same plot
plt.figure(figsize=(10, 6))

for i, autocorr in enumerate(autocorrelations):
    label = f'T ={temperatures[i]} K'
    plt.plot(time, autocorr, linestyle='-', label=label)

plt.axhline(y=1/np.e, color='red', linestyle='--', label='$1/e$')
plt.title('Autocorrelation Functions for Different Temperatures')
plt.xlabel('Time (t)')
plt.ylabel('$\chi(t)$')
plt.grid(True)
plt.legend()
plt.show()
