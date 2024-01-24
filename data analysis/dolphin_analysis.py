import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 

arguments = ['/Energy.txt', '/Magnetization.txt'] 
titles = [ 'Energy', 'Magnetization']

temperatures = [ 0.009 , 0.01 , 0.02 , 0.03 , 0.05 , 0.07 , 0.08 , 0.09 , 0.1 , 0.11 , 0.13 , 0.14 , 0.15, 0.2 , 0.25 , 0.3 , 0.35 , 
                0.4 , 0.45 , 0.5 , 0.55 , 0.6 , 0.65 , 0.7 , 0.75 , 0.8 , 0.85 , 0.9 , 0.92 , 0.94 , 0.95 , 0.96 , 0.98 , 1.0 , 1.02 , 1.04 , 1.05 , 1.07 ,
                  1.09 , 1.1 , 1.12 , 1.14 , 1.15 , 1.17 , 1.18 , 1.2 , 1.23 , 1.25 , 1.27 , 1.29 , 
                1.3 , 1.33 , 1.35 , 1.37 , 1.4 , 1.45 , 1.5 , 1.6 , 1.7 , 1.8 , 1.9 , 2.0 , 2.1 , 2.2 , 2.3 , 2.4 , 2.5 , 
                3.0 , 3.5 , 4.0 , 4.5 , 5.0  ]  

                

mean_argument_values = [] 
mean_helicity_values = []
mean_vortices_values = []

mean_vertex_values = []
mean_antivertex_values = []

N = 16 * 16



# Create a dictionary to store the arrays
mean_argument_values = {argument: [] for argument in arguments}

def calculate_std (data) : 
    try:

        std_deviation = np.std(data)
        return std_deviation
    
    except Exception as e:
        # Handle any exceptions that may occur during the calculation
        print(f"Error: {e}")
        return None


def calculate_mean (data) : 

    mean = sum(data) / len(data)
    return mean


print(len(temperatures)) 

i = 0
# Load the magnetization data and calculate autocorrelations for each temperature
for argument in arguments: 
    mean_argument_values[argument] = []

    for temp in temperatures  :
          # Create an empty list for each argument

        file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output/T_{temp}" + argument

        with open(file_path, 'r') as file:
     
            numbers = [float(line.strip()) for line in file.readlines()]
            #print(len(numbers))
            mm = calculate_mean (numbers)
            #print(mm)
            mean_argument_values[argument].append(mm)

    print(len(mean_argument_values[argument]))   

    temperature_matrix = np.column_stack((temperatures, mean_argument_values[argument]))

    plt.figure()  # Create a new figure for each plot
    plt.plot(temperature_matrix[:, 0],  temperature_matrix[:, 1], marker='o', linestyle='-')  
    plt.xlabel('Temperature (K)')
    plt.ylabel(titles[i])
    plt.title(f'Mean {titles[i]} vs. Temperatures (L=16, J = 1)')
    plt.grid(True)
    i += 1

plt.show()

#vortices and antivortices 

for temp in temperatures  :

    file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output/T_{temp}" + '/Vortices.txt'
    
    with open(file_path, 'r') as file:
     
        numbers = [float(line.strip()) for line in file.readlines()]
        mm = calculate_mean (numbers)
        #print(mm)
        mean_vertex_values.append(mm)
            
vertex_matrix = np.column_stack((temperatures, mean_vertex_values))

for temp in temperatures  :

    file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output/T_{temp}" + '/Antivortices.txt'

    with open(file_path, 'r') as file:
     
        numbers = [float(line.strip()) for line in file.readlines()]
        mm = calculate_mean (numbers)
        #print(mm)
        mean_antivertex_values.append(-mm)
            
antivertex_matrix = np.column_stack((temperatures, mean_antivertex_values))



plt.plot(vertex_matrix[:, 0], vertex_matrix[:, 1], marker='o',  linestyle='-', label= 'Vortex') #marker='o', 
plt.plot(antivertex_matrix[:, 0], antivertex_matrix[:, 1], marker='o', linestyle='-', label = 'Antivortex')
plt.legend()
plt.xlabel('Temperature (K)')
plt.ylabel('Vortex, Antivortex')
plt.title('Vortices and Antivortices vs. Temperatures (L=16, J = 1)')
plt.grid(True)

plt.show()


#Superfluid Stiffness

# Initialize empty arrays to store the columns
Jp = []
Js = []
mean_Jd_values = []

# Open the file for reading

i=0 

for temp in temperatures:

    Jd = []
    Ic = []

    file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output/T_{temp}" + '/Helicity_modulus.txt' 
    with open(file_path, "r") as file:

        # Iterate over each line in the file
        for line in file:
            
            columns = line.strip().split()
            if len(columns) == 2:
                
                Jd.append(float(columns[0]))
                Ic.append(float(columns[1]))

    mean_sin = calculate_std (Ic) ** 2
    
    sin = mean_sin * N / temperatures[i]

    Jp.append(sin)

    #print(sin)

    #squares = [num ** 2 for num in Ic]
    #mm_sq = calculate_mean (squares) / temperatures[i]

    mm = calculate_mean (Jd) 
    cos_Jd = mm 
    mean_Jd_values.append(cos_Jd)

    Js_new = cos_Jd - sin
    Js.append(Js_new)
    
    i += 1 
    #print(i)
    
Jd_matrix = np.column_stack((temperatures, mean_Jd_values))
Jp_matrix = np.column_stack((temperatures, Jp))
Js_matrix = np.column_stack((temperatures, Js))

y = [0] * len(temperatures)  # Initialize y with zeros

for temp in range(len(temperatures)): 
    y [temp] = 2 * temperatures [temp] / np.pi



plt.plot(Jd_matrix[:, 0], Jd_matrix[:, 1], linestyle='-', label = r'$J_d$') #marker='o', 
plt.plot(Jp_matrix[:, 0], Jp_matrix[:, 1],  linestyle='-', label = r'$J_p$')
plt.plot(Js_matrix[:, 0], Js_matrix[:, 1], linestyle='-', label = r'$J_s$')
plt.plot(temperatures , y, linestyle='--', label = r'$\frac{2T}{\pi}$')
plt.legend()
#plt.xlim(1, 2)
#plt.ylim(0, 0.2)
plt.xlabel('Temperature (K)')
plt.ylabel('superfluid stiffness')
plt.title('Superfulitity vs. Temperatures (L=16, J = 1)')
plt.grid(True)

plt.show()




