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
parser.add_argument('--L0', type=float, help='L0' )

# Parse the command-line arguments
args = parser.parse_args()

# Access the argument values
L = args.L
K = args.K
e = args.e
beta_high = args.b_high
beta_low = args.b_low
rank = args.rank
L0 = args.L0

# Now you can use these values in your script
print("L=", L)
print("K=", K)
print("e=", e)
print("beta high=", beta_high)
print("beta low=", beta_low)
print("rank=", rank)
print("L0=", L0)

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


#let us obtain the temperatures 
temperatures = []
delta = (1/beta_low - 1/beta_high)/(rank)
T_high = 1/beta_low
T_low = 1/beta_high
#print(delta)

i = 0

for l in L:
    
    N= l * l

    Jp1 = []
    Js1 = []
    mean_Jd_values1 = [] 
    sin_Ic1= []

    Js_inf = []

    Jp2 = []
    Js2 = []
    mean_Jd_values2 = [] 
    sin_Ic2= []

    helicity_sum_inf = []
    helicity = []
    temperatures = []

    for n in range(rank) : 

        t = T_high - n * delta
        print (t)
        
        temperatures.append(t)
        Jd1 = []
        Ic1 = []

        Jd2 = []
        Ic2 = []

        file_path1 = f"/home/x_mirpi/Output_TBG/K_{K}/e_{e}/L{l}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}/beta_{n}" + '/Helicity_modulus1.txt'
        with open(file_path1, "r") as file:
            for line in file:
                columns = line.strip().split()
                if len(columns) == 2:
                    Jd1.append(float(columns[0]))
                    Ic1.append(float(columns[1]))

        Jd1 = [x * N for x in Jd1]
        Ic1 = [x * N for x in Ic1]

        #Ic_std1 = calculate_std(Ic1) ** 2
        Ic_std1 = np.std(Ic1) ** 2
        sin_Ic1.append(Ic_std1)
        sin = Ic_std1 / t
        Jp1.append(sin)

        mm = calculate_mean(Jd1) 
        cos_Jd = mm 
        mean_Jd_values1.append(cos_Jd)

        Js_new1 = cos_Jd - sin
        Js1.append(Js_new1)

        file_path2 = f"/home/x_mirpi/Output_TBG/K_{K}/e_{e}/L{l}_K{K}_e{e}_bmin{beta_low}_bmax{beta_high}/beta_{n}" + '/Helicity_modulus2.txt'
        with open(file_path2, "r") as file:
            for line in file:
                columns = line.strip().split()
                if len(columns) == 2:
                    Jd2.append(float(columns[0]))
                    Ic2.append(float(columns[1]))

        Jd2 = [x * N for x in Jd2]
        Ic2 = [x * N for x in Ic2]

        #Ic_std2 = calculate_std(Ic2) ** 2
        Ic_std2 = np.std(Ic2) ** 2
        sin_Ic2.append(Ic_std2)
        #sin = Ic_std2 * N / t
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
        molt_mean = calculate_mean(Ic1)*calculate_mean(Ic2)
        mean_m_array = np.array(mean_molt)
        molt_m_array = np.array(molt_mean)
        sub = 1/t * (mean_m_array - molt_m_array)
        helicity_sum = (Js_new1 + Js_new2 -2*sub)/N
        #print(helicity_sum)

        helicity.append(helicity_sum)

        helicity_sum_crit = helicity_sum/ (1 + (1/(2 * np.log(l/L0))))
        helicity_sum_inf.append(helicity_sum_crit)

    Jd_matrix1 = np.column_stack((temperatures, mean_Jd_values1))
    Jp_matrix1 = np.column_stack((temperatures, Jp1))
    Js_matrix1 = np.column_stack((temperatures, helicity_sum_inf))
    
    plt.plot(Js_matrix1[:, 0], Js_matrix1[:, 1], label=f'L={l}', color=colors[i])

    i = i+1




y = [0] * len(temperatures) 

for temp in range(len(temperatures)): 
    y[temp] = 2 * temperatures[temp] / np.pi


plt.plot(temperatures , y, linestyle='--', label = r'$\frac{2T}{\pi}$', color = 'black')

plt.xlabel('T (K)')
plt.ylabel(r'$J_s$')
plt.title(f' $K={K}$ e = {e} $L_0 = {L0}$')
plt.legend()
plt.grid(True)
plt.savefig(f'Js_inf_K{K}_e{e}_L0{L0}.jpg')

plt.show()