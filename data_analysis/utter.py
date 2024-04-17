import numpy as np
from statsmodels.tsa.stattools import acf
import matplotlib.pyplot as plt

# Load the magnetization data from a file
magnetization_data = np.loadtxt("/Users/mirimi/Desktop/hihi/KTH/XY-model/Output/T_1.4/Magnetization.txt")
total_steps = len(magnetization_data)

# Calculate the autocorrelation using statsmodels
autocorrelation = acf(magnetization_data, nlags=250)
nlags=250

# Create lags for the x-axis
lags = np.arange(len(autocorrelation))
time = np.arange(len(autocorrelation)) * (total_steps / nlags)


# Plot the autocorrelation function
plt.figure(figsize=(10, 6))
plt.plot(lags, autocorrelation, marker='o', linestyle='-')
plt.title('Autocorrelation Function')
plt.xlabel('Lag')
plt.ylabel('Autocorrelation')
plt.grid(True)

# Threshold for estimating tau (e.g., 1/e)
threshold = 1 / 2.71828

# Find the lag at which autocorrelation drops below the threshold
tau = None
for lag, acf_value in enumerate(autocorrelation):
    if acf_value < threshold:
        tau = lag
        break
tau_time = tau * (total_steps / nlags)

if tau is not None:
    plt.axvline(x=tau_time, color='r', linestyle='--', label=f'Estimated τ = {tau}')
    plt.legend()

plt.show()
