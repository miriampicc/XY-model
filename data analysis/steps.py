import matplotlib.pyplot as plt
import numpy as np
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
parser.add_argument('--start', type=int, help='start')
parser.add_argument('--end', type=int, help='end')


args = parser.parse_args()

# Access the argument values
L = args.L
K = args.K
e = args.e
beta_high = args.b_high
beta_low = args.b_low
rank = args.rank
a = args.a
start = args.start
end = args.end

# Now you can use these values in your script
print("L=", L)
print("K=", K)
print("e=", e)
print("beta high=", beta_high)
print("beta low=", beta_low)
print("rank=", rank)
print("a=", a)
print("start=", start)
print("end=", end)

temperatures = []
delta = (1/beta_low - 1/beta_high)/(rank)
T_high = 1/beta_low
T_low = 1/beta_high

plt.figure(figsize=(10, 6))

beta_subset = 5

for n in range(start, end+1): 

    file_path = f"/home/x_mirpi/Output_TBG/K_{K}_tdf2/e_{e}/L{L}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}_a{a}/beta_{n}" + '/Rank.txt'

    with open(file_path, 'r') as file:
        data = [float(line.strip()) for line in file.readlines()]

    data = data[:20]

    # Generate the Monte Carlo step numbers (from 1 to 500000)
    monte_carlo_steps = np.arange(1, len(data) + 1)

    # Plot the step function

    plt.step(monte_carlo_steps, data, where='post', label=f'beta = {n}')

plt.xlabel('Monte Carlo Step')
plt.ylabel('Value')
plt.title('Monte Carlo Simulation Step Function')
plt.legend()

# Save the image to show as output
plt.savefig(f'Steps_a={a}.jpg')
#plt.show()