import matplotlib.pyplot as plt
import os

# 1. Definició de les dades (valors de la segona imatge)
x = [0.9, 1, 1.1, 1.2, 1.25, 1.275]

# Valors variables (punts blaus)
y_variable = [73.0887509, 73.7598943, 69.3147181, 64.1853886, 63.5988767, 63.5988767]

# 2. Creació de la gràfica
plt.figure(figsize=(10, 6))

# Pintem els punts blaus (variables) del Gas de Lennard-Jones
plt.scatter(x, y_variable, color='#1f77b4', label='Gas de Lennard-Jones')

# MODIFICACIÓ: Línia horitzontal constant i infinita per al Gas ideal 🚀
plt.axhline(y=69.3147181, color="#ff0e0e", linestyle='--', linewidth=2, label='Gas ideal')

# 3. Configuració d'eixos i títols
plt.xlabel('$T$ normalitzada')
plt.ylabel('$k_T$ normalitzada')
plt.grid(True, linestyle='-', alpha=0.3)
plt.legend(loc='upper right')

# Ajustem els límits adaptats a la nova imatge
plt.ylim(62, 76)
plt.xlim(0.85, 1.3)

# 4. MODIFICAT: Desa el gràfic a la carpeta 'Latex' 📁📈
directori_actual = os.path.dirname(os.path.abspath(__file__)) # .../Treball de Simulació/LennardJones
directori_pare = os.path.dirname(directori_actual)            # .../Treball de Simulació

# Ajuntem la ruta cap a la carpeta Latex
ruta_final = os.path.join(directori_pare, 'Latex', 'kt_densitat_gas.png')
plt.savefig(ruta_final)

# Opcional: mostrar-lo per pantalla
print(f"Gràfic generat i desat correctament a: {ruta_final}")
plt.show()