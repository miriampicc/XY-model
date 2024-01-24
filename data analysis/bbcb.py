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



for temp in temperatures : 
    file_path = f"/Users/mirimi/Desktop/hihi/KTH/XY-model/Output/T_{temp}" + '/Helicity_modulus.txt' 

    mm = 0

    Ic = []
    Jd = []

    with open(file_path, 'r') as file:
     
        for line in file:
            
                columns = line.strip().split()
                if len(columns) == 2:
                
                    Jd.append(float(columns[0]))
                    Ic.append(float(columns[1]))

    #print (Jd)
    mm = np.mean (Jd)
    print(f"Temperature: {temp}", mm)
    