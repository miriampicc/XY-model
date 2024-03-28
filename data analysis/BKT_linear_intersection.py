import numpy as np 
import pandas as pd
import math
import matplotlib.pyplot as plt 

numbers = ['1', '2']
arguments = ['/Magnetization1.txt', '/Magnetization2.txt'] 
title = ['Lattice 1', 'Lattice 2']


"""temperatures = [ 0.15, 
                0.17, 0.2, 0.22, 0.25, 0.27, 0.3, 0.32, 0.35, 0.37, 0.4, 0.45, 0.5, 0.52, 0.55, 
                0.57, 0.6, 0.62, 0.65, 0.67, 0.7, 0.72, 0.75, 0.77, 0.8,]  """

"""temperatures = [ 0.009, 0.01, 0.02, 0.03, 0.05, 0.07, 0.08, 0.09, 0.1, 0.11, 0.13, 0.14, 0.15, 
                0.17, 0.2, 0.22, 0.25, 0.27, 0.3, 0.32, 0.35, 0.37, 0.4, 0.45, 0.5, 0.52, 0.55, 
                0.57, 0.6, 0.62, 0.65, 0.67, 0.7, 0.72, 0.75, 0.77, 0.8, 0.82, 0.85, 0.87, 0.9, 
                0.92, 0.94, 0.95, 0.96, 0.98, 1.0, 1.02, 1.04, 1.05, 1.06, 1.07, 1.09, 1.1, 1.13, 1.15, 1.17, 1.2, 1.22, 1.25, 
                1.27, 1.3, 1.4, 1.45, 1.5, 1.7, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5]  """

"""temperatures =[0.3, 0.32, 0.35, 0.36, 0.37, 0.39, 0.4, 0.405, 0.41, 0.415, 0.42, 0.422, 0.424, 0.426, 0.428, 
               0.43, 0.432, 0.434, 0.436, 0.438, 0.44, 0.442, 0.444, 0.446, 0.448, 0.45, 0.455, 0.46, 0.465, 
               0.47, 0.475, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6,
                0.62, 0.65, 0.67, 0.7, 0.72, 0.75, 0.77, 0.8, 0.82, 0.85, 0.87, 0.9, 0.92, 0.94, 0.95, 0.96, 0.98, 1.0]"""

"""temperatures = [0.4, 0.405, 0.41, 0.415, 0.42, 0.422, 0.424, 0.426, 0.428, 0.43, 0.4305, 0.4310, 0.4315, 0.4320, 0.4325,
                0.4330, 0.4335, 0.4340, 0.4345, 0.4350, 0.4355, 0.4360, 0.4365, 0.4370, 0.4375, 0.4380, 0.4385, 0.4390, 0.4395, 0.44,
                0.4405, 0.4410, 0.4415, 0.4420, 0.4425, 0.4430, 0.4435, 0.4440, 0.4445, 0.4450, 0.4455, 0.4460, 0.4465, 0.4470, 0.4475,
                0.4480, 0.4485, 0.4490, 0.4495, 0.45, 0.4505, 0.4510, 0.4515, 0.4520, 0.4525, 0.4530, 0.4535, 0.4540, 0.4545, 0.4550,
                0.4555, 0.4560, 0.4565, 0.4570, 0.4575, 0.4580, 0.4585, 0.4590, 0.4595, 0.46]
"""

temperatures = [0.45, 0.455, 0.46, 0.465, 0.47, 0.475, 0.48, 0.485, 0.49, 0.495,
                0.5, 0.505, 0.51, 0.515, 0.52, 0.525, 0.53, 0.535, 0.54, 0.545,
                0.55, 0.555, 0.56, 0.565, 0.57, 0.575, 0.58, 0.585, 0.59, 0.595,
                0.6, 0.605, 0.61, 0.615, 0.62, 0.625, 0.63, 0.635, 0.64, 0.645,0.65]

"""temperatures = [0.5, 0.502, 0.504, 0.506, 0.508, 0.51, 0.512, 0.514, 0.516, 0.518, 0.52, 
                0.522, 0.524, 0.526, 0.528, 0.53, 0.532, 0.534, 0.536, 0.538, 0.54, 0.542, 
                0.544, 0.546, 0.548, 0.55, 0.552, 0.554, 0.556, 0.558, 0.56, 0.562, 0.564, 
                0.566, 0.568, 0.57, 0.572, 0.574, 0.576, 0.578, 0.58, 
                0.582, 0.584, 0.586, 0.588, 0.59, 0.592, 0.594, 0.596, 0.598, 0.6]"""

#L=[8, 12, 16, 20, 24, 32, 40, 48]  #, 20, 24, 32, 40
#L=[8, 12, 16, 20, 24, 32, 40, 48]
L=[8, 12, 16, 20, 24, 32, 40, 48]  #, 
rainbow_colors = ['red', 'orange', 'deeppink', 'indigo', 'violet', 'green', 'blue', 'brown' ]
K=5

# Create a dictionary to store the arrays
mean_argument_values = []

def calculate_std(data): 
    try:
        std_deviation = np.std(data)
        return std_deviation
    except Exception as e:
        # Handle any exceptions that may occur during the calculation
        print(f"Error: {e}")
        return None

def calculate_mean(data): 
    mean = sum(data) / len(data)
    return mean
b = 0 

for l in L: 
    # Superfluid Stiffness
    Jp1 = []
    Js1 = []
    mean_Jd_values1 = [] 
    sin_Ic1= []

    Jp2 = []
    Js2 = []
    mean_Jd_values2 = [] 
    sin_Ic2= []

    helicity_sum = []
    x_L = []
    temp_xL = []
    f = []

    N = l * l

    i = 0

    for temp in temperatures:
        Jd1 = []
        Ic1 = []

        Jd2 = []
        Ic2 = []

        #file_path1 = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Relative_fluctuations/K={K}/Output_L={l}/T_{temp}" + '/Helicity_modulus1.txt'
        file_path1 = f"/Users/mirimi/Desktop/hihi/KTH/XY-model-fixed/Output_K={K}/Output_L={l}/T_{temp}" + '/Helicity_modulus1.txt'
        with open(file_path1, "r") as file:
            for line in file:
                columns = line.strip().split()
                if len(columns) == 2:
                    Jd1.append(float(columns[0]))
                    Ic1.append(float(columns[1]))

        #file_path2 = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Relative_fluctuations/K={K}/Output_L={l}/T_{temp}" + '/Helicity_modulus2.txt'
        file_path2 = f"/Users/mirimi/Desktop/hihi/KTH/XY-model-fixed/Output_K={K}/Output_L={l}/T_{temp}" + '/Helicity_modulus2.txt'
        with open(file_path2, "r") as file:
            for line in file:
                columns = line.strip().split()
                if len(columns) == 2:
                    Jd2.append(float(columns[0]))
                    Ic2.append(float(columns[1]))

        mean_sin = calculate_std(Ic1) ** 2
        sin_Ic1.append(mean_sin)
        sin = mean_sin * N / temp

        Jp1.append(sin)

        mm = calculate_mean(Jd1) 
        cos_Jd = mm 
        mean_Jd_values1.append(cos_Jd)

        Js_new1 = cos_Jd - sin
        Js1.append(Js_new1)

        mean_sin = calculate_std(Ic2) ** 2
        sin_Ic2.append(mean_sin)
        sin = mean_sin * N / temp

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
        #print(mean_molt)
        molt_mean = calculate_mean(Ic1)*calculate_mean(Ic2)
        mean_m_array = np.array(mean_molt)
        molt_m_array = np.array(molt_mean)
        sub = 1/temp * (mean_m_array - molt_m_array)
        sott = Js_new1 + Js_new2 +2*sub
        #print (i, sub) 

        helicity_sum.append(sott) #important! To be put in the formula !!!!!!
        x_L_val = np.pi * (sott /(2*temp))
        print(x_L_val)
        if (x_L_val > 1 ): 
            x_L.append (x_L_val) 
            temp_xL.append(temp)
            f_val = math.log (l) - 1/(2*(x_L_val-1))
            print(f_val)
            f.append (f_val)

    #plt.plot(temperatures, helicity_sum, label = f'L={l}', color=rainbow_colors[b % len(rainbow_colors)])
    #plt.plot(temperatures, f, label = f'L={l}', color=rainbow_colors[b % len(rainbow_colors)])
    #plt.subplot(1,2,1)
    #plt.plot(temperatures, helicity_sum, label = f'L={l}', color=rainbow_colors[b % len(rainbow_colors)])
    #plt.xlabel('Temperature (K)')
    #plt.ylabel('Helicity Modulus sum')
    #plt.legend()
    #plt.grid(True)

    #plt.subplot(1,2,2)
    plt.scatter(temp_xL, f, color=rainbow_colors[b % len(rainbow_colors)])
    plt.plot(temp_xL, f, label = f'L={l}', color=rainbow_colors[b % len(rainbow_colors)])
    plt.xlabel('Temperature (K)')
    plt.ylabel(r'$f_L(T)=\ln(L)-\frac{1}{2}\left[\frac{1}{x_L(T)-1}\right]$')
    plt.legend()
    plt.grid(True)
    #plt.show()
    b += 1 

print(len(x_L))
print(len(temp_xL))
print(len(f))


plt.show()