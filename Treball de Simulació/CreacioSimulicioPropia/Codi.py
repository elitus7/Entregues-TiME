import random
import math
import matplotlib.pyplot as plt
import numpy as np

''' Apartat Q6 '''

N = 1000 # Nombre de partícules del sistema.
T = 200 # Temperatura (en K).

kb = 1.38e-23 # Constant de Boltzmann.
eps = 1e-21 # Valor d'epsilon que apareix a l'exercici 35 (en J).
beta = 1/(kb*T)

Nivells = [0,1*eps,10*eps] # Estats d'energia permesos per a cada partícula. 
Valors_E1 = []
Valors_E2 = []
Valors_E3 = []


PN = list(random.choice(Nivells) for _ in range(N)) # Posició de cada partícula a l'inici de la simulació (probabilitat de cada nivell de 1/3).

for i in range(200000): 
    i = random.randint(0,N-1)  # Selecciona de manera aleatoria una partícula.
    E_in = PN[i] # Definim E_in com l'energia inicial de la partícula escollida.
    E_nov = random.choice(Nivells) # Definim E_nov com l'energia en el nou nivell, que és escollit de manera aleatoria també.

    Dif_E = E_nov - E_in # Diferència d'energia en el canvi de nivell de la partícula.

    # Definim la condició següent (regla de Metropolis): Si aquesta diferència d'energia és inferior a 0, s'accepta el canvi de nivell. Si és superior a 0, es genera un número aleatori de 0,0 a 1,0. Si aquest és inferior al exp(-beta*Dif_E), s'accepta el canvi de nivell. Si no compleix cap de les condicions anteriors, la partícula es queda en el nivell inicial. 
    if Dif_E <= 0: 
        PN[i] = E_nov
    else:
        if random.random() < math.exp(-beta*Dif_E):
            PN[i] = E_nov

    # Cada 3000 passos registrem el nombre de partícules de cada nivell (per seprat), així més endevant podem comparar valors i confirmar que s'ha arribat a l'estat d'equilibri tèrmic. 
    if i % 3000 == 0: 
        Valors_E1.append(PN.count(0)) 
        Valors_E2.append(PN.count(1*eps))
        Valors_E3.append(PN.count(10*eps))

# Repetim el procés descrit - vegades per arribar a l'estat d'equilibri tèrmic. 


plt.hist(PN, bins=30, color='blue')
plt.xlim(-0.5 * eps, 11 * eps)
plt.show()


''' Apartat Q7 '''

Z = math.exp(-beta*0) + math.exp(-beta*eps) + math.exp(-beta*10*eps) # Funció de partició 

# AIXÒ SEGURAMENT HO BORRARÉ PERQUÈ HO HE ACABAT FENT MÉS GENERALMENT. 
# Per contestar l'apartat de l'estudi de l'ocupació mitajana de cada un dels tres nivells. Aquest és un cas particular (per a la T fixada al principi del codi). 

# Nombre de partícules al nivell 1 (E1 = 0 J), nivell 2 (E2 = 1*eps J) i nivell 3 (E3 = 10*eps J) respectivament que genera el codi. 
NP_Nivell1 = PN.count(0)
NP_Nivell2 = PN.count(1*eps)
NP_Nivell3 = PN.count(10*eps)

# Nombre de partícules al nivell 1 (E1 = 0 J), nivell 2 (E2 = 1*eps J) i nivell 3 (E3 = 10*eps J) respectivament que ens diu la teoria en l'estat d'equilibri tèrmic. 
NPT_Nivell1 = N*math.exp(-beta*0)/Z
NPT_Nivell2 = N*math.exp(-beta*eps)/Z
NPT_Nivell3 = N*math.exp(-beta*10*eps)/Z

# Comparació dels dos valors per cada nivell més l'error relatiu. 
print('Nivell 1: Numero de particules Codi= ', NP_Nivell1, 'Ocupacio teorica =', NPT_Nivell1, 'Error(%) =', abs(NP_Nivell1 - NPT_Nivell1)/NPT_Nivell1 )
print('NIvell 2: Numero de particules = ', NP_Nivell2, 'Ocupacio teorica =', NPT_Nivell2, 'Error(%) =', abs(NP_Nivell2 - NPT_Nivell2)/NPT_Nivell2 )
print('Nivell 3: Numero de particules = ', NP_Nivell3, 'Ocupacio teorica =', NPT_Nivell3, 'Error(%) =', abs(NP_Nivell3 - NPT_Nivell3)/NPT_Nivell3 )


''' Apartat Q7. Garantir que ens trobem en una situació d'equilibri. '''

# Dividim els valors obtinguts cada 3000 iteracions en dos. La primera part (2/3) la deixem com està. La segona part (1/3) li apliquem una regressió lineal i trobem el pendent d'aquesta. Ho fem per als tres nivells d'energia.
tall_1 = round((2/3*len(Valors_E1)))
tall_2 = round((2/3*len(Valors_E2)))
tall_3 = round((2/3*len(Valors_E3)))

eix_x = [1000*j for j in range(len(Valors_E1))] # Definim eix de les x per crear el gràfic. Serà igual per als tres nivells ja que tenim el mateix nombre de valors a les tres llistes Valors_Ei. 

# Definim els valors de x i y per als tres nivells per a la primera part (2/3).
x1_NE = eix_x[:tall_1]
y1_NE = Valors_E1[:tall_1]
x2_NE = eix_x[:tall_2]
y2_NE = Valors_E2[:tall_2]
x3_NE = eix_x[:tall_3]
y3_NE = Valors_E3[:tall_3]

# Definim els valors de x i y per als tres nivells per a la segona part (1/3).
x1_E = eix_x[tall_1:]
y1_E = Valors_E1[tall_1:]
x2_E = eix_x[tall_2:]
y2_E = Valors_E2[tall_2:]
x3_E = eix_x[tall_3:]
y3_E = Valors_E3[tall_3:]

# Regressions lineals per als tres nivells.
# Unitats dels pendents de [m] = Nombre de partícules/(3000*iteracions). Per arribar a l'equilibri, idealment [m] = 0. 
m1, n1 = np.polyfit(x1_E, y1_E, 1)
recta1 = [m1 * x + n1 for x in x1_E]
print('Pendent recta 1 = ', m1)

m2, n2 = np.polyfit(x2_E, y2_E, 1)
recta2 = [m2 * x + n2 for x in x1_E]
print('Pendent recta 2 = ',m2)

m3, n3 = np.polyfit(x3_E, y3_E, 1)
recta3 = [m3 * x + n3 for x in x3_E]
print('Pendent recta 3 = ', m3)

# Considerem equilibri tèrmic quan el pendent de l'ocupació (m) és inferior al llindar de tolerància (epsilon = 0.0005). És necessari que els tres pendents compleixin la condició. 
if abs(m1) < 0.0005 and abs(m2) < 0.0005 and abs(m3) < 0.0005:
    print("S'ha arribat a l'equilibri termic")

''' HE DE POSAR TÍTOLS ALS GRÀFICSSSSSSSS'''

# Gràfic Ocupació nivell 1 respecte cada iteració. A l'últim 1/3 de les dades s'aplica una regressió lineal. 
plt.plot(x1_NE, y1_NE, '.', color='blue')
plt.plot(x1_E, recta1, color='black', linewidth=2, label = 'Regressió lineal')
plt.plot(x1_E, y1_E, '.', color='red')
plt.ylim(min(Valors_E1) -20 , max(Valors_E1) + 20)
plt.title("Gràfic d'ocupació al Nivell 1 (0 J) respecte cada iteració ")
plt.xlabel("Iteracions")
plt.ylabel("Ocupació (partícules)")
plt.legend()
plt.show()

# Gràfic Ocupació nivell 2 respecte cada iteració. a l'últim 1/3 de les dades s'aplica una regressió lineal. 
plt.plot(x2_NE, y2_NE, '.', color='blue')
plt.plot(x2_E, recta2, color='black', linewidth=2, label = 'Regressió lineal')
plt.plot(x2_E, y2_E, '.', color='red')
plt.ylim(min(Valors_E2) -20 , max(Valors_E2) + 20)
plt.title("Gràfic d'ocupació al Nivell 2 (1*eps J) respecte cada iteració ")
plt.xlabel("Iteracions")
plt.ylabel("Ocupació (partícules)")
plt.legend()
plt.show()

# Gràfic Ocupació nivell 3 respecte cada iteració. a l'últim 1/3 de les dades s'aplica una regressió lineal. 
plt.plot(x3_NE, y3_NE, '.', color='blue')
plt.plot(x3_E, recta3, color='black', linewidth=2, label = 'Regressió lineal')
plt.plot(x3_E, y3_E, '.', color='red')
plt.ylim(min(Valors_E3) -20 , max(Valors_E3) + 20)
plt.title("Gràfic d'ocupació al Nivell 1 (10*eps J) respecte cada iteració ")
plt.xlabel("Iteracions")
plt.ylabel("Ocupació (partícules)")
plt.legend()
plt.show()

''' Q7. Estudi de l'ocupació mitjana de cada un dels tres nivells d'energia en l'equilibri en funció de la temperatura.'''

# Usant la funció de OcupacioCodi, els valors donats NO seran els mateixos que els obtinguts fora de la funció en cada ejecutació, tot i ser el mateix codi (degut a l'aleatorietat). Nomes la usem per a estudiar l'ocupació en funcio de la Temperatura i comparació amb el valor teòric. 

# Funció que ens retorna l'ocupació de cada nivell per a una certa temperatura basat en la distribució de Boltzmann.
def OcupacioTeorica(T): 
    beta = 1/(T*kb)
    Z = math.exp(-beta*0) + math.exp(-beta*eps) + math.exp(-beta*10*eps)

    OcupacioT = [N*math.exp(-beta*0)/Z, N*math.exp(-beta*eps)/Z, N*math.exp(-beta*10*eps)/Z]

    return OcupacioT

# Funció que ens retorna l'ocupació de cada nivell per a una certa temperatura generat per la regla de Metropolis (és exactament el codi del principi). 
def OcupacioCodi(T):
    beta = 1/(kb*T)

    PN = list(random.choice(Nivells) for _ in range(N)) # Posició de cada partícula a l'inici de la simulació (probabilitat de cada nivell de 1/3).

    for i in range(200000): 
        i = random.randint(0,N-1)  # Selecciona de manera aleatoria una partícula.
        E_in = PN[i] # Definim E_in com l'energia inicial de la partícula escollida.
        E_nov = random.choice(Nivells) # Definim E_nov com l'energia en el nou nivell, que és escollit de manera aleatoria també.

        Dif_E = E_nov - E_in # Diferència d'energia en el canvi de nivell de la partícula.

    # Definim la condició següent: Si aquesta diferència d'energia és inferior a 0, s'accepta el canvi de nivell. Si és superior a 0, es genera un número aleatori de 0,0 a 1,0. Si aquest és inferior al exp(-beta*Dif_E), s'accepta el canvi de nivell. Si no compleix cap de les condicions anteriors, la partícula es queda en el nivell inicial. 
        if Dif_E <= 0: 
            PN[i] = E_nov
        else:
            if random.random() < math.exp(-beta*Dif_E):
                PN[i] = E_nov
    
    OcupacioC = [PN.count(0), PN.count(1*eps), PN.count(10*eps)]

    return OcupacioC    

# DEFINICIÓ EIXOS X (Temperatura (K))

Valors_T = [20*x for x in range(1, 51)] # Valors de temperatura de 1 K a 1000 K amb intervals de 20 K (per a la simulació).
Valors_TT = np.linspace(1, 1000, 500) # Valors de temperatura de 1 K a 1000 K (per al càlcul analític).

# DEFINICIÓ EIXOS Y (Ocupacions (partícules))

# VALORS SIMULACIÓ
# Creem tres llistes buides per a calcular les ocupacions de cada nivell per a diferents temperatures (Valors_T).
Ocupacio_C_1 = [] 
Ocupacio_C_2 = [] 
Ocupacio_C_3 = [] 

for w in Valors_T:
    resultat_simulacio = OcupacioCodi(w)
    Ocupacio_C_1.append(resultat_simulacio[0]) # Ocupació en funció de la temperatura per al nivell 1 
    Ocupacio_C_2.append(resultat_simulacio[1]) # Ocupació en funció de la temperatura per al nivell 2
    Ocupacio_C_3.append(resultat_simulacio[2]) # Ocupació en funció de la temperatura per al nivell 3

# VALORS ANALÍTICS 
Ocupacio_T_1 = [OcupacioTeorica(x)[0] for x in Valors_TT] # Ocupació en funció de la temperatura per al nivell 1
Ocupacio_T_2 = [OcupacioTeorica(x)[1] for x in Valors_TT] # Ocupació en funció de la temperatura per al nivell 2
Ocupacio_T_3 = [OcupacioTeorica(x)[2] for x in Valors_TT] # Ocupació en funció de la temperatura per al nivell 3


# Gràfic comparatiu entre els resultats analítics (línia interlineada) i els resultats de la simulació (línia contínua) per als tres nivells d'energia. 

plt.plot(Valors_T, Ocupacio_C_1, color='blue', label = "Resultat Simulació, Nivell 1")
plt.plot(Valors_T, Ocupacio_C_2, color='red', label = "Resultat Simulació, Nivell 2")
plt.plot(Valors_T, Ocupacio_C_3, color='black', label = "Resultat Simulació, Nivell 3")

plt.plot(Valors_TT, Ocupacio_T_1 ,'--', color='blue', label = "Resultat Analític, Nivell 1")
plt.plot(Valors_TT, Ocupacio_T_2, '--', color='red', label = "Resultat Analític, Nivell 2")
plt.plot(Valors_TT, Ocupacio_T_3, '--', color='black', label = "Resultat Analític, Nivell 3")

plt.title("Comparativa resultats simulació amb resultats analítics.")
plt.ylabel("Ocupació (partícules)")
plt.xlabel("Temperatura (K)")
plt.legend()

plt.show()