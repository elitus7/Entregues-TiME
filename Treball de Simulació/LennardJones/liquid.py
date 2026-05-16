import matplotlib.pyplot as plt
import os

# 1. Definició de les dades de la nova imatge
x = [0.900, 1.000, 1.100, 1.200, 1.250, 1.275]
y_variable = [0.0789496, 0.09841019, 0.12198392, 0.12799133, 0.13181753, 0.13371857]

# 2. Creació de la gràfica
plt.figure(figsize=(10, 6))

# Pintem els punts blaus (variables)
plt.scatter(x, y_variable, color='#1f77b4', label='Líquid de Lennard-Jones')

# 3. Configuració d'eixos i títols
plt.xlabel('$T$ normalitzada (adimensional)')
plt.ylabel('$k_T$ normalitzada (adimensional)')
plt.grid(True, linestyle='-', alpha=0.3)
plt.legend(loc='lower right') # Canviat a la dreta-baix perquè no tapi els punts

# Ajustem els límits de la vista clònics a la teva imatge de Excel
plt.ylim(0.07, 0.14)
plt.xlim(0.85, 1.3)

# 4. Desa el gràfic a la carpeta 'Latex'
directori_actual = os.path.dirname(os.path.abspath(__file__)) 
directori_pare = os.path.dirname(directori_actual)            

# Guardem el fitxer amb un nom representatiu de la temperatura
ruta_final = os.path.join(directori_pare, 'Latex', 'kt_liquid.png')
plt.savefig(ruta_final)

# Mostrar per pantalla
print(f"Gràfic 'kt_liquid.png' generat i desat correctament a: {ruta_final}")
plt.show()