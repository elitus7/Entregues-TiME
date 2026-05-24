import math
import matplotlib.pyplot as plt
import numpy as np

# Constants i variables. Per aquestes últimes, posem valors arbitraris però raonables. 
eps_kb = 175.12  # epsilon/kB      
rho = 1.7023              
r0_cm = 3.6396e-8          
Na = 6.02214e23           

A = (2*np.pi*Na*(r0_cm**3))/3  

T_aprox = eps_kb/(0.06) # Temperatura llindar a la que comença a ser vàlida l'aproximació (Error inferior al 3%).
 
# Funció B_2(T) exacta (trobada a l'apartat b).
def fB2(T):
    beta = 1/T
    Resultat = A*(rho**3 + math.exp(beta*eps_kb)*(1-rho**3))

    return Resultat

# Funció B_2(T) aproximada amb sèrie de Taylor. 
def fB2_aprox(T): 
    beta = 1/T
    Resultat = A*(1 + beta*eps_kb*(1-rho**3))

    return Resultat

# EIX X (Temperatura (K))
Temp = [150 + 2*x for x in range(1, 2000)]

# EIX Y (B_2(T) exacta i aproximada (cm^3/mol))
B2 = []
B2_aprox = []

for i in Temp:
    valors = fB2(i)
    B2.append(valors)

for j in Temp: 
    valors = fB2_aprox(j)
    B2_aprox.append(valors)


# Gràfica comparativa entre B_2(T) exacta i aproximada. Indicant la temperatura llindar i les zones vàlides i no vàlides per l'aproximació. 
plt.plot(Temp, B2, '--', color='red', label='$B_2(T)$ Exacte')
plt.plot(Temp, B2_aprox, color='black', label='$B_2(T)$ Aproximat')
plt.axvline(x=T_aprox, color='green', linestyle=':', linewidth=2, label='$T$ llindar')
plt.axvspan(min(Temp), T_aprox, color='red', alpha=0.07, label='Zona de no validesa (Error > 3%)')
plt.axvspan(T_aprox , max(Temp), color='green', alpha=0.07, label='Zona de validesa (Error $\leq$ 3%)')
plt.xlabel("Temperatura (K)")
plt.ylabel(r'$B_2$ ($\text{cm}^3/\text{mol}$)')
plt.ylim(min(B2_aprox), max(B2_aprox) + 20)
plt.xlim(min(Temp), max(Temp)) # Ajusta l'eix X exactament al teu rang
plt.grid(True, alpha=0.3)
plt.legend(loc='lower right') 
plt.show()

# Temperatura de Boyle i d'inversió màxima calculades a partir de B_2(T) aproximada
T_B = eps_kb*(rho**3 - 1)
T_inv = 2*eps_kb*(rho**3 - 1)

# Gràfica de B_2(T) aproximada més T_B i T_inv
plt.plot(Temp, B2_aprox, color='black', label='$B_2(T)$ Aproximat')
plt.axvline(x=T_B, color='green', linestyle=':', linewidth=2, label='$T$ de Boyle aproximada')
plt.axvline(x=T_inv, color='red', linestyle=':', linewidth=2, label="$T$ d'inversió màxima aproximada")
plt.axvline(x=T_aprox, color='blue', linestyle=':', linewidth=2, label='$T$ llindar')
plt.xlabel("Temperatura (K)")
plt.ylabel(r'$B_2$ ($\text{cm}^3/\text{mol}$)')
plt.xlim(0 , max(Temp))
plt.ylim(min(B2_aprox), max(B2_aprox) + 20)
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
