import numpy as np
import matplotlib.pyplot as plt

# Define your magnetization values M(t)
# Replace with your own data

magnetization_data = np.loadtxt("/Users/mirimi/Desktop/hihi/KTH/XY-model/Output/T_2.2/Magnetization.txt") 
M_t = magnetization_data[:1000]

avg_magn = np.mean(magnetization_data)

# Perform FFT to obtain M(ω)
M_omega = np.fft.fft(M_t)

M_omega[0] = 0

# Calculate |M(ω)|^2 to obtain χ(ω)
chi_omega = np.abs(M_omega) ** 2

# Inverse FFT to obtain χ(t)
chi_t = np.fft.ifft(chi_omega)

# Normalize the functions
M_omega_norm = np.abs(M_omega) / np.max(np.abs(M_omega))
chi_omega_norm = chi_omega / np.max(chi_omega)
chi_t_norm = chi_t / np.max(chi_t)

# Calculate the corresponding frequencies
N = len(M_t)  # Number of data points
frequencies = np.fft.fftfreq(N)

# Plot the results for visualization
plt.figure(figsize=(12, 4))

# Plot M(ω)
plt.subplot(131)
plt.plot(frequencies, M_omega_norm)
plt.xlabel('Frequency (Hz)')
plt.ylabel('|M(ω)|')
plt.title('M(ω)')
plt.grid()

# Plot χ(ω)
plt.subplot(132)
plt.plot(frequencies, chi_omega_norm)
plt.xlabel('Frequency (Hz)')
plt.ylabel('χ(ω)')
plt.title('χ(ω)')
plt.grid()

# Plot χ(t)
plt.subplot(133)
plt.plot(np.arange(0, N), chi_t_norm)
plt.xlabel('Time')
plt.ylabel('χ(t)')
plt.title('χ(t)')
plt.grid()

plt.tight_layout()
plt.show()
