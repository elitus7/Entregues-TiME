import matplotlib.pyplot as plt
import os

# 1. Definició de les dades
x = [0.9, 1, 1.1, 1.2, 1.25, 1.275]

# Valors de dalt (variables)
y_variable = [90.0847414, 83.4623758, 82.5201305, 82.3083279, 81.5869273, 81.5735469]

# 2. Creació de la gràfica
plt.figure(figsize=(10, 6))

# Pintem els punts blaus (variables)
plt.scatter(x, y_variable, color='#1f77b4', label='Gas de Lennard-Jones')

# SOLUCIÓ: Dibuixa la línia constant de punta a punta absoluta del gràfic 🚀
plt.axhline(y=81.0930216, color="#fa0000", linestyle='--', linewidth=2, label='Gas ideal')

# 3. Configuració d'eixos i títols
plt.xlabel('$T$ normalitzada (adimensional)')
plt.ylabel('$k_T$ normalitzada (adimensional)')
plt.grid(True, linestyle='-', alpha=0.3)
plt.legend(loc='upper right')

# Ajustem els límits de la vista
plt.ylim(80, 92)
plt.xlim(0.85, 1.3)

# 4. MODIFICAT: Desa el gràfic a la carpeta 'Latex' 📁📈
directori_actual = os.path.dirname(os.path.abspath(__file__)) # .../Treball de Simulació/LennardJones
directori_pare = os.path.dirname(directori_actual)            # .../Treball de Simulació

# Ajuntem la ruta cap a la carpeta Latex
ruta_final = os.path.join(directori_pare, 'Latex', 'kt_volum_gas.png')
plt.savefig(ruta_final)

# Opcional: mostrar-lo per pantalla
print(f"Gràfic 'kt_volum_gas.png' generat i desat correctament a: {ruta_final}")
plt.show()