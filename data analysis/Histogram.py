from collections import Counter
import matplotlib.colors as mcolors
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 
import argparse

parser = argparse.ArgumentParser(description='Description of the script')

# Add arguments
parser.add_argument('--L', nargs='+', type=int, help='L')
parser.add_argument('--K', type=float, help='K')
parser.add_argument('--e', type=float, help='e')
parser.add_argument('--b_high', type=float, help='beta high')
parser.add_argument('--b_low', type=float, help='beta low')
parser.add_argument('--rank', type=int, help='rank')
parser.add_argument('--beta', type=int, help='beta')


# Parse the command-line arguments
args = parser.parse_args()

# Access the argument values
L = args.L
K = args.K
e = args.e
beta_high = args.b_high
beta_low = args.b_low
rank = args.rank
beta = args.beta

# Now you can use these values in your script
print("L=", L)
print("K=", K)
print("e=", e)
print("beta high=", beta_high)
print("beta low=", beta_low)
print("rank=", rank)
print("beta =", beta)



def calculate_std(data): 
    try:
        std_deviation = np.std(data)
        std_deviation = std_deviation / (np.sqrt(len(data)-1))
        return std_deviation
    except Exception as e:
        # Handle any exceptions that may occur during the calculation
        print(f"Error: {e}")
        return None

def calculate_mean(data): 
    mean = sum(data) / len(data)
    return mean

def plot_histogram(data):
    # Count frequency of each number from 0 to 63
    counter = Counter(data)
    frequencies = [counter[i] for i in range(64)]

    # Plotting

    start_color = 'red'
    end_color = 'lightblue'
    cmap = mcolors.LinearSegmentedColormap.from_list('red_to_blue', [start_color, end_color], N=64)


    colors = [cmap(i / 63) for i in range(64)]

    fig, ax = plt.subplots()

    bars = ax.bar(range(64), frequencies, tick_label=range(64), color=colors)
    ax.set_xlabel('Rank')
    ax.set_ylabel('Frequency')
    ax.set_title(f'Frequency of moving of the replica {beta}')

    # Set x-ticks to multiples of 5
    ax.set_xticks(range(0, 64, 5))

    # Create a ScalarMappable for the colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=63))
    sm.set_array([])

    # Add the colorbar to the plot, stealing space from the axes
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Rank')


    plt.savefig(f'Histogram_e={e}_K={K}_beta{beta}.jpg')
    plt.show()


hist = []
betas = []
delta_beta = (beta_high - beta_low)/(rank)

for l in L: 
    N = l * l

    #file_path = f"/Users/mirimi/Desktop/OUTPUT_cluster/e_{e}/L{l}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}/beta_{n}" + '/Energy.txt'
    file_path = f"/home/x_mirpi/Output_TBG/K_{K}/e_{e}/L{l}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}/beta_{beta}" + '/Rank.txt'
    #print(n)

    with open(file_path, 'r') as file:
        hist = [float(line.strip()) for line in file.readlines()]

    plot_histogram(hist)

