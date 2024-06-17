from collections import Counter
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

colors = [
    "#00008B",  # Dark Blue
    "#0000CD",  # Medium Blue
    "#4169E1",  # Royal Blue
    "#1E90FF",  # Dodger Blue
    "#00BFFF",  # Deep Sky Blue
    "#87CEEB",  # Sky Blue
    "#87CEFA",  # Light Sky Blue
    "#EEE8AA",  # Pale Goldenrod
    "#F0E68C",  # Khaki
    "#FFD700",  # Gold
    "#F08080",  # Light Coral
    "#FF6347"   # Tomato
]

def plot_histogram(data):
    # Count frequency of each number from 0 to 63
    counter = Counter(data)
    frequencies = [counter[i] for i in range(64)]

    # Plotting
    plt.bar(range(64), frequencies, tick_label=range(64))
    plt.xlabel('Rank')
    plt.ylabel('Frequency')
    plt.title('Frequency of Numbers from 0 to 63')
    plt.savefig(f'Histogram_e={e}_K={K}.jpg')
    plt.show()


hist = []
betas = []
delta_beta = (beta_high - beta_low)/(rank)

for l in L: 
    N = l * l

    #file_path = f"/Users/mirimi/Desktop/OUTPUT_cluster/e_{e}/L{l}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}/beta_{n}" + '/Energy.txt'
    file_path = f"/home/x_mirpi/Output_TBG/K_{K}/e_{e}/L{L}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}/beta_{beta}" + '/Energy.txt'
    #print(n)

    with open(file_path, 'r') as file:
        hist = [float(line.strip()) for line in file.readlines()]

    plot_histogram(hist)

