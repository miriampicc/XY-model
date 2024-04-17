import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 
from shapely.geometry import LineString
from scipy.optimize import curve_fit

numbers = ['1', '2']
arguments = ['/Magnetization1.txt', '/Magnetization2.txt'] 
title = ['Lattice 1', 'Lattice 2']


"""temperatures = [ 0.009, 0.01, 0.02, 0.03, 0.05, 0.07, 0.08, 0.09, 0.1, 0.11, 0.13, 0.14, 0.15, 
                0.17, 0.2, 0.22, 0.25, 0.27, 0.3, 0.32, 0.35, 0.37, 0.4, 0.45, 0.5, 0.52, 0.55, 
                0.57, 0.6, 0.62, 0.65, 0.67, 0.7, 0.72, 0.75, 0.77, 0.8, 0.82, 0.85, 0.87, 0.9, 
                0.92, 0.94, 0.95, 0.96, 0.98, 1.0, 1.02, 1.04, 1.05, 1.06, 1.07, 1.09, 1.1, 1.13, 1.15, 1.17, 1.2, 1.22, 1.25, 
                1.27, 1.3, 1.4, 1.45, 1.5] """  
"""
temperatures =[0.009, 0.01, 0.02, 0.03, 0.05, 0.07, 0.08 ,0.09, 0.1, 0.11, 0.13, 0.14, 0.15, 0.17, 0.2, 0.22, 0.25, 
                0.27, 0.29, 0.3, 0.32, 0.35, 0.36, 0.37, 0.39, 0.4, 0.42, 0.44, 0.45, 0.46, 0.47, 0.49, 0.5, 0.51, 
                0.52, 0.53, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6]"""

"""temperatures =[0.3, 0.32, 0.35, 0.36, 0.37, 0.39, 0.4, 0.405, 0.41, 0.415, 0.42, 0.422, 0.424, 0.426, 0.428, 
               0.43, 0.432, 0.434, 0.436, 0.438, 0.44, 0.442, 0.444, 0.446, 0.448, 0.45, 0.455, 0.46, 0.465, 
               0.47, 0.475, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6,
                0.62, 0.65, 0.67, 0.7, 0.72, 0.75, 0.77, 0.8, 0.82, 0.85, 0.87, 0.9, 0.92, 0.94, 0.95, 0.96, 0.98, 1.0]"""

temperatures = [0.4, 0.405, 0.41, 0.415, 0.42, 0.422, 0.424, 0.426, 0.428, 0.43, 0.4305, 0.4310, 0.4315, 0.4320, 0.4325,
                0.4330, 0.4335, 0.4340, 0.4345, 0.4350, 0.4355, 0.4360, 0.4365, 0.4370, 0.4375, 0.4380, 0.4385, 0.4390, 0.4395, 0.44,
                0.4405, 0.4410, 0.4415, 0.4420, 0.4425, 0.4430, 0.4435, 0.4440, 0.4445, 0.4450, 0.4455, 0.4460, 0.4465, 0.4470, 0.4475,
                0.4480, 0.4485, 0.4490, 0.4495, 0.45, 0.4505, 0.4510, 0.4515, 0.4520, 0.4525, 0.4530, 0.4535, 0.4540, 0.4545, 0.4550,
                0.4555, 0.4560, 0.4565, 0.4570, 0.4575, 0.4580, 0.4585, 0.4590, 0.4595, 0.46]


L = [8, 12, 20]

#L=[8, 16, 32, 64]
#L=[8, 12]
rainbow_colors = ['red', 'orange', 'deeppink', 'green', 'blue', 'indigo', 'violet', 'brown']
K=1
intersections = []
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

def logarithmic_function(x, a, b):
    return a * np.log(x) + b

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

    double = []

    N = l * l

    i = 0

    for temp in temperatures:
        Jd1 = []
        Ic1 = []

        Jd2 = []
        Ic2 = []

        file_path1 = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Relative_fluctuations/K={K}_bis/Output_L={l}/T_{temp}" + '/Helicity_modulus1.txt'
        with open(file_path1, "r") as file:
            for line in file:
                columns = line.strip().split()
                if len(columns) == 2:
                    Jd1.append(float(columns[0]))
                    Ic1.append(float(columns[1]))

        file_path2 = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Relative_fluctuations/K={K}_bis/Output_L={l}/T_{temp}" + '/Helicity_modulus2.txt'
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

        double.append(sott)
    

    
    b += 1 


    y = [0] * len(temperatures)

    for temp in range(len(temperatures)): 
        y[temp] = 2 * temperatures[temp] / np.pi


    Line1 = LineString (np.column_stack((temperatures, double)))
    Line2 = LineString (np.column_stack((temperatures, y)))
    intersection = Line1.intersection (Line2)
    intersections.append(intersection)

    plt.plot(temperatures, double, label = f'L={l}', color=rainbow_colors[b % len(rainbow_colors)])
    plt.plot(*intersection.xy, 'ro')

    print(intersection)

plt.plot(temperatures , y, linestyle='--', label = r'$\frac{2T}{\pi}$', color='black')
plt.legend()
plt.xlabel('Temperature (K)')
plt.ylabel('Helicity sum')
plt.title('Helicity sum of bilayer compound')
plt.grid(True)

# Add a title to the entire figure
#plt.suptitle('($J_1 = J_2 = K$)')

#plt.tight_layout()  # Adjust layout for better appearance



plt.show()

print (intersections)

# Separate x and y components of intersections
intersection_x = []
intersection_y = []

for intersection in intersections:
    intersection_x.append(intersection.x)
    intersection_y.append(intersection.y)

int_x = np.array(intersection_x)
int_y = np.array(intersection_y) 

for i in range(len(int_x)): 

    plt.plot(intersection_x[i], intersection_y[i], 'o',label=f'L={L[i]}', color=rainbow_colors[i % len(rainbow_colors)])

plt.xlabel('T')
plt.ylabel(r'$J_s$')
plt.title('Intersection points')
plt.legend()
plt.grid()
plt.show()



#Function of the INVERSE lattice dimensions
b=0
plt.plot(1/L[0], int_x[0], 'o', label=f'L={L[0]}', color='red' )
plt.plot(1/L[1], int_x[1], 'o', label=f'L={L[1]}', color='orange' )
plt.plot(1/L[2], int_x[2], 'o', label=f'L={L[2]}', color='deepskyblue' )
plt.plot(1/L[3], int_x[3],'o', label=f'L={L[3]}', color='deeppink' )
plt.plot(1/L[4], int_x[4],'o', label=f'L={L[4]}', color='violet' )
plt.plot(1/L[5], int_x[5],'o', label=f'L={L[5]}', color='indigo' )
plt.plot(1/L[6], int_x[6],'o', label=f'L={L[6]}', color='brown' )

plt.grid()
plt.xlabel(r'$\frac{1}{L}$')
plt.ylabel(r'$T_C$')
plt.grid(True)
plt.legend()
plt.title('Critical T as function of the lattice dim')
plt.show()

#linear fitting 
x_data = np.array([1/L[0], 1/L[1], 1/L[2], 1/L[3], 1/L[4], 1/L[5], 1/L[6]])   #
y_data = np.array([int_x[0], int_x[1], int_x[2], int_x[3], int_x[4], int_x[5], int_x[6]])  #

#x_data = np.array([1/L[0], 1/L[1], 1/L[2], 1/L[3]]) 
#y_data = np.array([int_x[0], int_x[1], int_x[2], int_x[3]]) 

fit_coefficients = np.polyfit(x_data, y_data, 1)
fit_function = np.poly1d(fit_coefficients)
print(fit_coefficients)
print(fit_function)

plt.plot(1/L[0], int_x[0], 'o', label=f'L={L[0]}', color='red' )
plt.plot(1/L[1], int_x[1], 'o', label=f'L={L[1]}', color='orange' )
plt.plot(1/L[2], int_x[2], 'o', label=f'L={L[2]}', color='deepskyblue' )
plt.plot(1/L[3], int_x[3],'o', label=f'L={L[3]}', color='green' )
plt.plot(1/L[4], int_x[4],'o', label=f'L={L[4]}', color='violet' )
plt.plot(1/L[5], int_x[5],'o', label=f'L={L[5]}', color='indigo' )
plt.plot(1/L[6], int_x[6],'o', label=f'L={L[6]}', color='brown' )
plt.plot(x_data, fit_function(x_data), linestyle='--',label='Linear Fit',linewidth=0.8, color='black')

plt.grid()
plt.xlabel(r'$\frac{1}{L}$')
plt.ylabel(r'$T_C^{BKT}$')
plt.legend()
plt.title('Critical T as a function of lattice dimension')
plt.show()


#function of the lattice dimensions
b=0

plt.plot(L[0], int_x[0], 'o', label=f'L={L[0]}', color='red' )
plt.plot(L[1], int_x[1], 'o', label=f'L={L[1]}', color='orange' )
plt.plot(L[2], int_x[2], 'o', label=f'L={L[2]}', color='deepskyblue' )
plt.plot(L[3], int_x[3],'o', label=f'L={L[3]}', color='green' )
plt.plot(L[4], int_x[4],'o', label=f'L={L[4]}', color='violet' )
plt.plot(L[5], int_x[5],'o', label=f'L={L[5]}', color='indigo' )
plt.plot(L[6], int_x[6],'o', label=f'L={L[6]}', color='brown' )

plt.grid()
plt.xlabel(r'$L$')
plt.ylabel(r'$T_C^{BKT}$')
plt.grid(True)
plt.legend()
plt.title('Critical T as function of the lattice dim')
plt.show()




