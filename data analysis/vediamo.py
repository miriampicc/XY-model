import numpy as np
from statsmodels.tsa.stattools import acf
import matplotlib.pyplot as plt

# Load the magnetization data from a file
magnetization_data = np.loadtxt("/Users/mirimi/Desktop/hihi/KTH/XY-model/Output/T_1.4/Magnetization.txt")

# Calculate the autocorrelation using statsmodels
autocorrelation = acf(magnetization_data, nlags=50)

# Create lags for the x-axis
lags = np.arange(len(autocorrelation))

# Plot the autocorrelation function
plt.figure(figsize=(10, 6))
plt.plot(lags, autocorrelation, marker='o', linestyle='-')
plt.title('Autocorrelation Function')
plt.xlabel('Lag')
plt.ylabel('Autocorrelation')
plt.grid(True)
plt.show()
