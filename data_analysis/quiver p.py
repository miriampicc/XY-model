import numpy as np
import matplotlib.pyplot as plt

argument = '/Spin_dir_fin.txt'

def angle_to_vector(angles):
    return np.array([np.cos(angles), np.sin(angles)])


def plot_spins(angles):
    vectors = angle_to_vector(angles)

    x, y = np.meshgrid(np.arange(grid_size[0]), np.arange(grid_size[1]))
    
    plt.quiver(x, y, vectors[0], vectors[1], angles, pivot='middle') #cmap='hsv'
    #plt.colorbar()
    plt.title('Lattice at T = 1.2 K (L=16, J = 1)')
    plt.show()

cos_sin = []

file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output/T_1.2" + argument 

with open(file_path, 'r') as file:
            angles = [float(line.strip()) for line in file.readlines()]

grid_size = (16, 16)
plt.grid(True)

plot_spins(angles)