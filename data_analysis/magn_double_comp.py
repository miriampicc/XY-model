import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt 

numbers = ['1', '2']
arguments = ['/Magnetization1.txt', '/Magnetization2.txt'] 

temperatures = [ 0.009 , 0.01 , 0.02 , 0.03 , 0.05 , 0.07 , 0.08 , 0.09 , 0.1 , 0.11 , 0.13 , 0.14 , 0.15, 0.17 , 0.2 , 0.22 , 0.25 , 0.27 , 0.3 , 0.32 , 0.35 , 0.37 ,
                0.4 , 0.45 , 0.5 , 0.55 , 0.6 , 0.65 , 0.7 , 0.75 , 0.85 , 0.9 , 0.95 , 1.05 ,
                  1.1 , 1.15 , 1.2 , 1.25 , 1.27 , 
                1.3 , 1.4 , 1.45 , 1.5 , 1.7 , 2.0 , 2.5 , 
                3.0 , 3.5 , 4.0 , 4.5 , 5.0  ]  

N = 16 * 16

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


# Create a dictionary to store the arrays
mean_argument_values = {argument: [] for argument in arguments}

#PLOT OF THE MAGNETISATION FIRST LAYER
i = 0
for argument in arguments: 

    mean_argument_values[argument] = []

    print(argument)

    for temp in temperatures:
        

        file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output/T_{temp}" + argument

        with open(file_path, 'r') as file:
            numbers = [float(line.strip()) for line in file.readlines()]
            mm = calculate_mean(numbers)
            mean_argument_values[argument].append(mm) 

        i += 1
        print(i)
   

for argument, values in mean_argument_values.items():
    plt.plot(temperatures, values, marker='o', linestyle='-', label=f'Argument: {argument}')

plt.xlabel('Temperature (K)')
plt.ylabel('Magnetisation')
plt.title('Magnetisation of the two layers')
plt.grid(True)
plt.legend()
plt.show()