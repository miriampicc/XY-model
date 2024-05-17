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

# Parse the command-line arguments
args = parser.parse_args()

# Access the argument values
L = args.L
K = args.K
e = args.e
beta_high = args.b_high
beta_low = args.b_low
rank = args.rank

# Now you can use these values in your script
print("L=", L)
print("K=", K)
print("e=", e)
print("beta high=", beta_high)
print("beta low=", beta_low)
print("rank=", rank)

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

#let us obtain the temperatures 
temperatures = []
delta = (1/beta_low - 1/beta_high)/(rank)
T_high = 1/beta_low
T_low = 1/beta_high
#print(delta)

for l in L: 

    i = 0

    mean_argument_values = {0: [], 1: []}
    std_values = {0: [], 1: []}
    temperatures = []


    N = l * l
 
    for n in range(rank):

        t = T_high - n * delta
        print (t)
        
        temperatures.append(t)

        layer1_values = []
        layer2_values = []
        mm1=0
        mm2=0
        std1=0
        std2=0

        #file_path = f"/home/x_mirpi/Output_TBG/K_{K}_first/e_{e}/L{l}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}/beta_{n}" + '/Single_Magnetization.txt'
        file_path = f"/Users/mirimi/Desktop/OUTPUT_cluster/e_{e}/L{l}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}/beta_{n}" + '/Single_Magnetization.txt'


        with open(file_path, 'r') as file:

            for line in file:
                    values = [float(val) for val in line.strip().split()]
                    layer1_values.append(values[0])
                    layer2_values.append(values[1])

            print(len(layer1_values))

            mm1 = calculate_mean(layer1_values)
            mm2 = calculate_mean(layer2_values)
            std1 = calculate_std (layer1_values)
            std2 = calculate_std (layer2_values)

            mean_argument_values[0].append(mm1)
            mean_argument_values[1].append(mm2)
            std_values[0].append(std1)
            std_values[1].append(std2)        
    
    magn1 = np.array(mean_argument_values[0])
    magn2 = np.array (mean_argument_values[1])
    m_std1 = np.array(std_values[0])
    m_std2 = np.array(std_values[1])

    plt.plot(temperatures, mean_argument_values[0], label=f'L={l}')   #, color=rainbow_colors[i % len(rainbow_colors)]
    plt.fill_between(temperatures, magn1 - m_std1, magn1 + m_std1, alpha=0.3, linewidth=4)

    i += 1

    plt.plot(temperatures, mean_argument_values[1])
    #plt.fill_between(temperatures, magn2 - m_std2, magn2 + m_std2, alpha=0.6, facecolor='lightblue', linewidth=4, label='Error Bounds ')


plt.xlabel('T')
plt.ylabel(r'$m$')
plt.title(f'Magnetization of single layers')
plt.legend()
plt.grid(True)
plt.savefig(f'Magnetization_single_layer_e={e}_K={K}.jpg')

plt.show()