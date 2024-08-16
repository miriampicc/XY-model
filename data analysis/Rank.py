import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import argparse

parser = argparse.ArgumentParser(description='Description of the script')

# Add arguments
parser.add_argument('--L', type=int, help='L')
parser.add_argument('--K', type=float, help='K')
parser.add_argument('--e', type=float, help='e')
parser.add_argument('--b_high', type=float, help='beta high')
parser.add_argument('--b_low', type=float, help='beta low')
parser.add_argument('--rank', type=int, help='rank')
parser.add_argument('--a', type=float, help='a')
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
a = args.a
beta = args.beta

# Now you can use these values in your script
print("L=", L)
print("K=", K)
print("e=", e)
print("beta high=", beta_high)
print("beta low=", beta_low)
print("rank=", rank)
print("a=", a)
print("beta=", beta)

delta = (1 / beta_low - 1 / beta_high) / rank
T_high = 1 / beta_low
T_low = 1 / beta_high
file_path = f'/home/x_mirpi/Output_TBG/K_{K}_tdf2/e_{e}/L{L}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}_a{a}/beta_{beta}/Rank.txt'

data = np.loadtxt(file_path, dtype=int)

# Step 2: Count the occurrences
counter = Counter(data)

# Sort by number (0 to 63) to ensure all bins are accounted for
numbers = list(range(64))
occurrences = [counter.get(number, 0) for number in numbers]

# Step 3: Plot the histogram
plt.bar(numbers, occurrences, color='blue')
plt.xlabel('Number')
plt.ylabel('Occurrences')
plt.xticks(range(0, 65, 5))
plt.title('Histogram of Number Occurrences')
plt.xticks(numbers)  # Ensure all numbers 0-63 are labeled
plt.grid(True)
plt.savefig(f'Rank_e={e}_replica={beta}.jpg')

plt.show()
