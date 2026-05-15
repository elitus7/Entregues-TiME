import random
import math
import matplotlib.pyplot as plt
import numpy as np


''' Apartat Q6. '''

N = 1000 # Nombre de partícules del sistema.
T = 10 # Temperatura (en K).
It = 200000 # Nombre d'iteracions. Si al executar salta un missatge 'Cal augmentar el nombre d'iteracions', cal augmentar aquest valor. 

kb = 1.38e-23 # Constant de Boltzmann.
eps = 1e-21 # Valor d'epsilon que apareix a l'exercici 35 (en J). Valor arbitrari. 
beta = 1/(kb*T)

Nivells = [0,1*eps,10*eps] # Estats d'energia permesos per a cada partícula. 
Valors_E1 = []
Valors_E2 = []
Valors_E3 = []


PN = list(random.choice(Nivells) for _ in range(N)) # Posició de cada partícula a l'inici de la simulació (probabilitat de cada nivell de 1/3).

for _ in range(It): 
    indexPart = random.randint(0,N-1)  # Selecciona de manera aleatoria una partícula.
    E_in = PN[indexPart] # Definim E_in com l'energia inicial de la partícula escollida.
    E_nov = random.choice(Nivells) # Definim E_nov com l'energia en el nou nivell, que és escollit de manera aleatoria també.

    Dif_E = E_nov - E_in # Diferència d'energia en el canvi de nivell de la partícula.

    # Definim la condició següent (regla de Metropolis): Si aquesta diferència d'energia és inferior a 0, s'accepta el canvi de nivell. Si és superior a 0, es genera un número aleatori de 0,0 a 1,0. Si aquest és inferior al exp(-beta*Dif_E), s'accepta el canvi de nivell. Si no compleix cap de les condicions anteriors, la partícula es queda en el nivell inicial. 
    if Dif_E <= 0: 
        PN[indexPart] = E_nov
    else:
        if random.random() <= math.exp(-beta*Dif_E):
            PN[indexPart] = E_nov

    # Cada uns certs passos registrem el nombre de partícules de cada nivell (per seprat), així més endevant podem comparar valors i confirmar que s'ha arribat a l'estat d'equilibri tèrmic. 
    if _ % round(It/100) == 0: 
        Valors_E1.append(PN.count(0)) 
        Valors_E2.append(PN.count(1*eps))
        Valors_E3.append(PN.count(10*eps))

# Repetim el procés descrit (It) vegades per arribar a l'estat d'equilibri tèrmic.


# Histograma d'ocupació a l'equilibri tèrmic.
Nhist= [PN.count(0),PN.count(1*eps), PN.count(10*eps)]
Etiq = [0,1,10]

'''

plt.bar(Etiq, Nhist, color='skyblue', edgecolor='black')
plt.xticks(Etiq) # Només mostra els números 0, 1 i 10
plt.xlabel(r"Nivell energètic ($\epsilon$ J)")
plt.title(f"Histograma d'ocupació a l'equilibri tèrmic a {T} K")
plt.ylabel("Ocupació (partícules)")
plt.show()
'''

''' Apartat Q7. Estudi de l'ocupació mitjana de cada un dels tres nivells d'energia en l'equilibri en funció de la temperatura.'''

# Atenció: Usant la funció de OcupacioCodi, els valors donats NO seran els mateixos que els obtinguts fora de la funció en cada executació, tot i ser el mateix codi (degut a l'aleatorietat). Nomes la usem per a estudiar l'ocupació en funcio de la Temperatura i comparació amb el valor teòric. 

# Funció que ens retorna l'ocupació de cada nivell per a una certa temperatura basat en la distribució de Boltzmann.
def OcupacioTeorica(T): 

    beta = 1/(T*kb)
    Z = math.exp(-beta*0) + math.exp(-beta*eps) + math.exp(-beta*10*eps) # Funció de partició 

    OcupacioT = [N*math.exp(-beta*0)/Z, N*math.exp(-beta*eps)/Z, N*math.exp(-beta*10*eps)/Z]

    return OcupacioT

# Funció que ens retorna l'ocupació de cada nivell per a una certa temperatura generat per la regla de Metropolis (S'ha adaptat el codi original integrant una mitjana temporal (100 mesures separades per 1000 passos) per minimitzar les fluctuacions i millorar la convergència visual de la simulació.). 
def OcupacioCodi(T):
    beta = 1/(kb*T)
    PN = [random.choice(Nivells) for _ in range(N)]
    
    # Primer bucle fa el mateix que el del codi principal. 
    for _ in range(It): 
        idx = random.randint(0, N-1)
        E_in = PN[idx]
        E_nov = random.choice(Nivells)
        if (E_nov - E_in) <= 0 or random.random() < math.exp(-beta*(E_nov - E_in)):
            PN[idx] = E_nov
            
    # Un cop arribat a l'equilibri, es mesura 100 cops l'estat a l'equilibri i se'n fa la mitjana dels valors d'ocupació per reduir el soroll. 
    m_n0, m_n1, m_n2 = 0, 0, 0
    for _ in range(100):
        for _ in range(1000):
            idx = random.randint(0, N-1)
            E_in = PN[idx]
            E_nov = random.choice(Nivells)
            if (E_nov - E_in) <= 0 or random.random() < math.exp(-beta*(E_nov - E_in)):
                PN[idx] = E_nov
        
        m_n0 += PN.count(0)
        m_n1 += PN.count(eps)
        m_n2 += PN.count(10*eps)
  
    return [m_n0/100, m_n1/100, m_n2/100] # Ens retorna l'ocupació mitjana a l'equilibri dels tres nivells d'energia. 

# DEFINICIÓ EIXOS X (Temperatura (K))
'''
Valors_T = [20*x for x in range(1, 151)] # Valors de temperatura de 1 K a 1000 K amb intervals de 20 K (per a la simulació).
Valors_TT = np.linspace(1, 3000, 500) # Valors de temperatura de 1 K a 1000 K (per al càlcul analític).

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
Ocupacio_T_3 = [OcupacioTeorica(x)[2]for x in Valors_TT] # Ocupació en funció de la temperatura per al nivell 3


# Gràfic comparatiu entre els resultats analítics (línia interlineada) i els resultats de la simulació (línia contínua) per als tres nivells d'energia. 

plt.figure(figsize=(10, 6))

plt.plot(Valors_T, Ocupacio_C_1, '.', color='blue', label = "Resultat Simulació, Nivell 1")
plt.plot(Valors_T, Ocupacio_C_2, '.',color='red', label = "Resultat Simulació, Nivell 2")
plt.plot(Valors_T, Ocupacio_C_3, '.',color='black', label = "Resultat Simulació, Nivell 3")

plt.plot(Valors_TT, Ocupacio_T_1 ,'--', color='blue', label = "Resultat Analític, Nivell 1")
plt.plot(Valors_TT, Ocupacio_T_2, '--', color='red', label = "Resultat Analític, Nivell 2")
plt.plot(Valors_TT, Ocupacio_T_3, '--', color='black', label = "Resultat Analític, Nivell 3")

plt.title("Comparativa resultats simulació amb resultats analítics.")
plt.ylabel("Ocupació (partícules)")
plt.xlabel("Temperatura (K)")
plt.legend()

plt.show()

'''
''' Per un cas particular (a T fixada) '''

# Aquest és un cas particular (per a la T fixada al principi del codi). Trobem l'error relatiu entre el resultat de la simulació i el resultat analític. Per a diferents resultats de T, només cal canviar el valor de l'inici i executar.

# Usem les funcions definides anteriorment. 
OcupC = OcupacioCodi(T)
OcupT = OcupacioTeorica(T)

# Nombre de partícules al nivell 1 (E1 = 0 J), nivell 2 (E2 = 1*eps J) i nivell 3 (E3 = 10*eps J) respectivament que genera la simulació.
NP_Niv1 = OcupC[0]
NP_Niv2 = OcupC[1]
NP_Niv3 = OcupC[2]

# Nombre de partícules al nivell 1 (E1 = 0 J), nivell 2 (E2 = 1*eps J) i nivell 3 (E3 = 10*eps J) respectivament basat en la distribució de Boltzmann.
NPT_Niv1 = OcupT[0]
NPT_Niv2 = OcupT[1]
NPT_Niv3 = OcupT[2]

# Comparació dels dos resultats per cada nivell més l'error relatiu. 
print('Nivell 1: Numero de particules Codi= ', NP_Niv1, 'Ocupacio teorica =', NPT_Niv1, 'Error(%) =', abs(NP_Niv1 - NPT_Niv1)*100/NPT_Niv1 )
print('NIvell 2: Numero de particules = ', NP_Niv2, 'Ocupacio teorica =', NPT_Niv2, 'Error(%) =', abs(NP_Niv2 - NPT_Niv2)*100/NPT_Niv2 )
print('Nivell 3: Numero de particules = ', NP_Niv3, 'Ocupacio teorica =', NPT_Niv3, 'Error(%) =', abs(NP_Niv3 - NPT_Niv3)*100/NPT_Niv3 )


''' Apartat Q7. Garantir que ens trobem en una situació d'equilibri (Cas particular, a T fixada). '''

# Per comprovar que s'ha arribat a l'estat d'equilibri tèrmic, previament s'ha guardat les dades de l'ocupació cada unes certes iteracions. 
# Fem una regressió lineal a una part final d'aquestes dades. Si el pendent d'aquesta regressió compleix el criteri de convergència establert, considerarem que s'ha arribat a l'equilibri (idealment, caldria que el pendent fos de 0, ja que les dades no haurien de variar a l'equilibri). 

# Dividim els valors obtinguts en dos. La primera part (2/3) la deixem com està. La segona part (1/3) fem una regressió lineal i trobem el pendent d'aquesta. Ho fem per als tres nivells d'energia.
tall_1 = round((2/3*len(Valors_E1))) # Fem que el tall es situi a 2/3 dels valors (per als tres nivells).
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

# Regressions lineals per als tres nivells. m1, m2 i m3 son els pendents de les regressions lineals que volem calcular per als tres nivells energètics respectivament. 
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
    print("S'ha arribat a l'equilibri tèrmic")

else: 
    print("Cal augmentar el nombre d'iteracions.")

# Escrivim les tres equacions de la recta per incloure-les a les gràfiques 
eq1 = f"y = {m1:.5f}x + {n1:.5f}"
eq2 = f"y = {m2:.5f}x + {n2:.5f}"
eq3 = f"y = {m3:.5f}x + {n3:.5f}"

# Gràfic Ocupació nivell 1 respecte cada iteració. A l'últim 1/3 de les dades s'aplica la regressió lineal. 
plt.plot(x1_NE, y1_NE, '.', color='blue')
plt.plot(x1_E, recta1, color='black', linewidth=2, label = 'Regressió lineal', zorder=10)
plt.plot(x1_E, y1_E, '.', color='red', zorder=5)
plt.text(max(eix_x)/2 , (max(Valors_E1)+min(Valors_E1))/2 , eq1, fontsize=9, color='black', 
         bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray'))
plt.ylim(min(Valors_E1) -20 , max(Valors_E1) + 20)
plt.title("Gràfic d'ocupació al Nivell 1 (0 J) respecte cada iteració ")
plt.xlabel("Iteracions")
plt.ylabel("Ocupació (partícules)")
plt.legend()
plt.show()

# Gràfic Ocupació nivell 2 respecte cada iteració. a l'últim 1/3 de les dades s'aplica la regressió lineal. 
plt.plot(x2_NE, y2_NE, '.', color='blue')
plt.plot(x2_E, recta2, color='black', linewidth=2, label = 'Regressió lineal', zorder=10)
plt.plot(x2_E, y2_E, '.', color='red', zorder=5)
plt.text(max(eix_x)/2 , min(Valors_E2) , eq2, fontsize=9, color='black', 
         bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray'))
plt.ylim(min(Valors_E2) -20 , max(Valors_E2) + 20)
plt.title("Gràfic d'ocupació al Nivell 2 (1*eps J) respecte cada iteració ")
plt.xlabel("Iteracions")
plt.ylabel("Ocupació (partícules)")
plt.legend()
plt.show()

# Gràfic Ocupació nivell 3 respecte cada iteració. a l'últim 1/3 de les dades s'aplica una regressió lineal. 
plt.plot(x3_NE, y3_NE, '.', color='blue')
plt.plot(x3_E, recta3, color='black', linewidth=2, label = 'Regressió lineal', zorder=10)
plt.plot(x3_E, y3_E, '.', color='red', zorder=5)
plt.text(max(eix_x)/2 , (max(Valors_E3)+min(Valors_E3))/2 , eq3, fontsize=9, zorder=10,
         bbox=dict(facecolor='white', alpha=0.8))
plt.ylim(min(Valors_E3) -20 , max(Valors_E3) + 20)
plt.title("Gràfic d'ocupació al Nivell 3 (10*eps J) respecte cada iteració ")
plt.xlabel("Iteracions")
plt.ylabel("Ocupació (partícules)")
plt.legend()
plt.show()

''' Q7. Valor de Tc '''

Tc = 10*eps/(kb*math.log(N)) # Temperatura crítica trobada al exercici 35. 
Tc_m = [] # Llista on s'afegiràn els valors de temperatura crítica trobats. S'ha optat per fer una mitjana de varis valors perquè la Tc trobava sigui exacta. 

# Per trobar Tc, s'ha optat per separar el càlcul del bloc principal per tal de disminuir l'interval entre els valors de temperatura, 
# així obtenint una major presició sense que el temps d'executació agumenti considerablement (ja que per estudiar l'ocupació en funció de la temperatura no ens fa falta un interval tan precís).

Valors_T2 = [0.5*x for x in range(120,300)] # Rang de temperatures més petit i més precís (+- 0,5) (ja que la Tc teòrica és d'uns 105) comparat amb el bloc principal (que va de 1 K a 3000 K(intervals de 20 K)).

for _ in range(10): # Trobem 10 vegades la temperatura crítica. Es farà una mitjana dels valors. 

    Tc_sim = None # Definim Tc que es calcularà amb la simulació.
    for k in Valors_T2:
        ResultatsSim = OcupacioCodi(k)

    if ResultatsSim[2] >= 1 and Tc_sim is None: 
        Tc_sim = k
        Tc_m.append(Tc_sim)
        break  # Un cop s'ha trobat la temperatura crítica, el bucle s'anul·la. 

print("Temperatura crítica simulació = ", sum(Tc_m)/10)
print("Temperatura crítica teòrica = ", Tc)
print("Error relatiu (%) = ", abs(sum(Tc_m)/10 - Tc)*100/Tc)

	

''' Apartat Q8. '''

# ESTÀ MALAMENT
'''
def Energia(N):
    PN = list(random.choice(Nivells) for _ in range(N)) # Posició de cada partícula a l'inici de la simulació (probabilitat de cada nivell de 1/3).

    for i in range(20000): 
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
    
    EnergiaEq = PN.count(0)*0 + PN.count(1*eps)*eps + PN.count(10*eps)*10*eps

    return EnergiaEq

Part = [10*x for x in range(1,500)]

Energy = [1/Energia(x) for x in Part]

plt.plot(Part, Energy, color='blue') 
plt.show()

'''
