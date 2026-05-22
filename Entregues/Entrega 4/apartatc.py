import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os
import matplotlib.ticker as ticker

# Dades experimentals de la taula (T en K, B2 en cm^3/mol)
T_data = np.array([273.16, 298.16, 323.16, 348.16, 373.16, 423.16, 473.16, 573.16, 673.16, 773.16, 873.16, 973.16])
B2_data = np.array([-154.74, -130.21, -110.62, -95.04, -82.13, -62.10, -46.74, -25.06, -10.77, -0.13, 7.95, 14.22])

# B2 expressió teòrica
def b2_teoric(T, eps_k, rho, r0_cm):
    N_A = 6.02214e23
    return (2 * np.pi * N_A * (r0_cm**3)) / 3 * (rho**3 + np.exp(eps_k / T) * (1 - rho**3))

# Valors inicials (posarem valors raonables perquè acabi convergint)
valors_inicials = [200.0, 1.5, 4e-8]

# Ajust de les dades experimentals a l'expressió teòrica
popt, pcov = curve_fit(b2_teoric, T_data, B2_data, p0=valors_inicials)
eps_k_opt, rho_opt, r0_opt = popt

# Visualitzaió dels resultats
print("Resultats de l'ajust:")
print(f"epsilon / k_B = {eps_k_opt:.2f} K")
print(f"rho (r1/r0)   = {rho_opt:.4f}")
print(f"r0            = {r0_opt:.4e} cm")

# Gràfica
plt.figure(figsize=(8, 6))
plt.scatter(T_data, B2_data, color='red', label='Dades experimentals')
T_fit = np.linspace(250, 1000, 200)
ax = plt.gca()
decimals = ticker.FormatStrFormatter('%.2f')
ax.xaxis.set_major_formatter(decimals)
ax.yaxis.set_major_formatter(decimals)
plt.plot(T_fit, b2_teoric(T_fit, *popt), label='Ajust teòric', color='blue')
plt.xlabel('Temperatura (K)')
plt.ylabel(r'$B_2$ ($cm^3/mol$)')
plt.legend()
plt.grid(True)

directori_actual = os.path.dirname(os.path.abspath(__file__))
directori= os.path.dirname(directori_actual)            

guardar = os.path.join(directori, 'Entrega 4', 'AjustB2.png')
plt.savefig(guardar, dpi=300, bbox_inches='tight')
plt.show()