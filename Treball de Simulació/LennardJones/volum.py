import matplotlib.pyplot as plt

# 1. Definició de les dades
x = [0.9, 1, 1.1, 1.2, 1.25, 1.275]

# Valors de dalt (variables)
y_variable = [90.0847414, 83.4623758, 82.5201305, 82.3083279, 81.5869273, 81.5735469]

# Valors de baix (constants)
y_constant = [81.0930216] * len(x)

# 2. Creació de la gràfica
plt.figure(figsize=(10, 6))

# Pintem els punts blaus (variables) i els taronja (constants)
plt.scatter(x, y_variable, color='#1f77b4', label='Gas de Lennard-Jones')
plt.scatter(x, y_constant, color="#fa0000", label='Gas ideal')

# 3. Configuració d'eixos i títols
plt.title('Compressibilitat isotèrmica $k_T$ en funció de la temperatura $T$')
plt.xlabel('$T$ normalitzada (adimensional)')
plt.ylabel('$k_T$ normalitzada (adimensional)')
plt.grid(True, linestyle='-', alpha=0.3)
plt.legend(loc='upper right')

# Ajustem els límits perquè es vegi semblant a la teva imatge
plt.ylim(80, 92)
plt.xlim(0.85, 1.3)

# 4. Desa el gràfic (això sobrescriurà el fitxer anterior cada cop que s'executi)
plt.savefig('kt_volum.png')

# Opcional: mostrar-lo per pantalla
print("Gràfic 'kt_volum.png' generat correctament!")
plt.show()