import numpy as np 

rank = 64 

T_low = 0.555
T_high = 0.625 

delta = (T_high - T_low)/ rank 

temps = []

for n in range(rank) : 

    t = T_low + n * delta 

    temps. append (t) 

    #print (t)

temps.reverse()

formatted_temps = " ".join(f"{temp:.3f}" for temp in temps)
print(formatted_temps)

print(len(temps))
    
