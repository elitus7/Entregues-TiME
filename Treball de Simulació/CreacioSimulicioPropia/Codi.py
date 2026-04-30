import random
import math
import matplotlib.pyplot as plt

N = 1000 # Nombre de partícules del sistema.
T = 10 # Temperatura (en K).

kb = 1.38e-23 # Constant de Boltzmann.
eps = 1e-21 # Valor d'epsilon que apareix a l'exercici 35 (en J).
beta = 1/(kb*T)

Nivells = [0,1*eps,10*eps] # Estats d'energia permesos per a cada partícula. 

PN = list(random.choice(Nivells) for _ in range(N)) # Posició de cada partícula a l'inici de la simulació (probabilitat de cada nivell de 1/3).


for i in range(7000): 
    i = random.randint(0,N-1)  # Selecciona de manera aleatoria una partícula.
    E_in = PN[i] # Definim E_in com l'energia inicial de la partícula escollida.
    E_nov = random.choice(Nivells) # Definim E_nov com l'energia en el nou nivell, que és escollit de manera aleatoria també.

    Dif_E = E_nov - E_in # Diferència d'energia en el canvi de nivell de la partícula.

    if Dif_E <= 0: # Definim la condició següent: Si aquesta diferència d'energia és inferior a 0, s'accepta el canvi de nivell. Si és superior a 0, es genera un número aleatori de 0,0 a 1,0. Si aquest és inferior al exp(-beta*Dif_E), s'accepta el canvi de nivell. Si no compleix cap de les condicions anteriors, la partícula es queda en el nivell inicial. 
        PN[i] = E_nov
    else:
        if random.random() < math.exp(-beta*Dif_E):
            PN[i] = E_nov

# Repetim el procés descrit 5000 vegades per arribar a l'estat d'equilibri tèrmic. 

plt.hist(PN, bins=30, color='blue')
plt.show()
