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


# Create a dictionary to store the arrays
mean_argument_values = []

def calculate_std(data): 
    try:

        std_deviation = np.mean(data**2) - (np.mean(data))**2
        #std_deviation = np.std(data)
        return std_deviation
    except Exception as e:
        # Handle any exceptions that may occur during the calculation
        print(f"Error: {e}")
        return None

def calculate_mean(data): 
    mean = sum(data) / len(data)
    return mean


#let us obtain the temperatures 

delta = (1/beta_low - 1/beta_high)/(rank)
T_high = 1/beta_low
T_low = 1/beta_high

# Superfluid Stiffness


for l in L :

    i = 0 
    Jp1 = []
    Js1 = []
    mean_Jd_values1 = [] 
    sin_Ic1= []

    Jp2 = []
    Js2 = []
    mean_Jd_values2 = [] 
    sin_Ic2= []

    double = []
    temperatures = []

    N = l * l 

    for n in range(rank):

        t = T_high - n * delta
        print (t)
            
        temperatures.append(t)

        Jd1 = []
        Ic1 = []

        Jd2 = []
        Ic2 = []

        #file_path1 = f"/home/x_mirpi/Output_TBG/K_{K}_tdf2/e_{e}/L{l}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}/beta_{n}" + '/Helicity_modulus1.txt'
        file_path1 = f"/home/x_mirpi/Output_TBG/K_{K}_tdf2/e_{e}/L{l}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}_a{a}/beta_{n}" + '/Helicity_modulus1.txt'


        with open(file_path1, "r") as file:
            for line in file:
                columns = line.strip().split()
                if len(columns) == 2:
                    Jd1.append(float(columns[0]))
                    Ic1.append(float(columns[1]))

        
        #file_path2 = f"/home/x_mirpi/Output_TBG/K_{K}_tdf2/e_{e}/L{l}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}/beta_{n}" + '/Helicity_modulus2.txt'
        file_path2 = f"/home/x_mirpi/Output_TBG/K_{K}_tdf2/e_{e}/L{l}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}_a{a}/beta_{n}" + '/Helicity_modulus2.txt'

        with open(file_path2, "r") as file:
            for line in file:
                columns = line.strip().split()
                if len(columns) == 2:
                    Jd2.append(float(columns[0]))
                    Ic2.append(float(columns[1]))

        Jd1 = [x * N for x in Jd1]
        Ic1 = [x * N for x in Ic1]

        Ic_std1 = np.std(Ic1) ** 2
        sin_Ic1.append(Ic_std1)
        sin = Ic_std1 / t
        Jp1.append(sin)


        mm = calculate_mean(Jd1)
        cos_Jd = mm
        mean_Jd_values1.append(cos_Jd)

        Js_new1 = cos_Jd - sin
        Js1.append(Js_new1)

        Jd2 = [x * N for x in Jd2]
        Ic2 = [x * N for x in Ic2]

        Ic_std2 = np.std(Ic2) ** 2
        sin_Ic2.append(Ic_std2)
        sin = Ic_std2 / t

        Jp2.append(sin)

        mm = calculate_mean(Jd2)
        cos_Jd = mm
        mean_Jd_values2.append(cos_Jd)

        Js_new2 = cos_Jd - sin
        Js2.append(Js_new2)

        result = []
        for i in range(len(Ic1)):
            result.append(Ic1[i] * Ic2[i])
        
        mean_molt = calculate_mean (result)
        print(mean_molt)
        molt_mean = calculate_mean(Ic1)*calculate_mean(Ic2)
        mean_m_array = np.array(mean_molt)
        molt_m_array = np.array(molt_mean)
        sub = 1/(t) * (mean_m_array - molt_m_array)
        sott = (Js_new1 + Js_new2 - 2*sub)/N
        #print (i, sub)

        double.append(sott)

        
        i += 1 

    Jd_matrix1 = np.column_stack((temperatures, mean_Jd_values1))
    Jp_matrix1 = np.column_stack((temperatures, Jp1))
    Js_matrix1 = np.column_stack((temperatures, Js1))

    Jd_matrix2 = np.column_stack((temperatures, mean_Jd_values2))
    Jp_matrix2 = np.column_stack((temperatures, Jp2))
    Js_matrix2 = np.column_stack((temperatures, Js2))

    duouble_matrix = np.column_stack ((temperatures, double))

y = [0] * len(temperatures)

for temp in range(len(temperatures)): 
    y[temp] = 2 * temperatures[temp] / np.pi


plt.figure(figsize=(15, 5))

# Plot 1
plt.subplot(1, 3, 1)  # 3 rows, 1 column, plot 1
plt.plot(Jd_matrix1[:, 0], Jd_matrix1[:, 1], label='Jd')
plt.plot(Jp_matrix1[:, 0], Jp_matrix1[:, 1], label='Jp')
plt.plot(Js_matrix1[:, 0], Js_matrix1[:, 1], label='Js')
plt.plot(temperatures , y, linestyle='--', label = r'$\frac{2T}{\pi}$')

plt.xlabel('Temperature (K)')
plt.ylabel('superfluid stiffness')
plt.title(' Component 1')
plt.legend()
plt.grid(True)

# Plot 2
plt.subplot(1, 3, 2)  # 3 rows, 1 column, plot 2
plt.plot(Jd_matrix2[:, 0], Jd_matrix2[:, 1], label='Jd')
plt.plot(Jp_matrix2[:, 0], Jp_matrix2[:, 1], label='Jp')
plt.plot(Js_matrix2[:, 0], Js_matrix2[:, 1], label='Js')
plt.plot(temperatures , y, linestyle='--', label = r'$\frac{2T}{\pi}$')

plt.legend()
plt.xlabel('Temperature (K)')
plt.ylabel('superfluid stiffness')
plt.title('Component 2')
plt.grid(True)

# Plot 3
plt.subplot(1, 3, 3)  # 3 rows, 1 column, plot 3
plt.plot(duouble_matrix[:, 0], duouble_matrix[:, 1])
plt.legend()
plt.xlabel('Temperature (K)')
plt.ylabel('Helicity sum')
plt.title('Helicity sum of bilayer compound')
plt.grid(True)

# Add a title to the entire figure
plt.suptitle(f'Total density fluctuations K={K}, e={e}, L={l}')

plt.tight_layout()  # Adjust layout for better appearance
plt.savefig(f'Multiple_superflidity_bla_bla_K{K}_e{e}.jpg')



plt.show()
