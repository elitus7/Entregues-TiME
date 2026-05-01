from vpython import *
import numpy as np
import math
import random

win = 500
Natoms = 500 # Nombre d'àtoms ajustat a 200 para rendimiento.

# Paràmetres de la simulació (tots els valors estan en SI).
L = 1 
gray = color.gray(0.7) 
mass = 4E-3/6E23 
Ratom = 0.04 
k = 1.4E-23 
T = 300 
dt = 1E-5

animation = canvas(width=win, height=win, align='left')
animation.range = L
animation.title = 'Hard Sphere Gas (Energy Distribution)'
s = """
Theoretical and averaged energy distributions (Joules).
Initially all atoms have the same energy, but collisions change them.
One of the atoms is marked and leaves a trail.
"""
animation.caption = s

d = L/2+Ratom
r = 0.005
boxbottom = curve(color=gray, radius=r)
boxbottom.append([vector(-d,-d,-d), vector(-d,-d,d), vector(d,-d,d), vector(d,-d,-d), vector(-d,-d,-d)])
boxtop = curve(color=gray, radius=r)
boxtop.append([vector(-d,d,-d), vector(-d,d,d), vector(d,d,d), vector(d,d,-d), vector(-d,d,-d)])
vert1 = curve(color=gray, radius=r)
vert2 = curve(color=gray, radius=r)
vert3 = curve(color=gray, radius=r)
vert4 = curve(color=gray, radius=r)
vert1.append([vector(-d,-d,-d), vector(-d,d,-d)])
vert2.append([vector(-d,-d,d), vector(-d,d,d)])
vert3.append([vector(d,-d,d), vector(d,d,d)])
vert4.append([vector(d,-d,-d), vector(d,d,-d)])

Atoms = []
p = [] 
apos = [] 
pavg = sqrt(2*mass*1.5*k*T) # Energia cinètica promig: p**2/(2mass) = (3/2)kT.

for i in range(Natoms):
    x = L*random.random()-L/2
    y = L*random.random()-L/2
    z = L*random.random()-L/2
    if i == 0:
        Atoms.append(sphere(pos=vector(x,y,z), radius=Ratom, color=color.cyan, make_trail=True, retain=100, trail_radius=0.3*Ratom))
    else:
        Atoms.append(sphere(pos=vector(x,y,z), radius=Ratom, color=gray))
    apos.append(vec(x,y,z)) 
    theta = pi*random.random() 
    phi = 2*pi*random.random()
    px = pavg*sin(theta)*cos(phi) 
    py = pavg*sin(theta)*sin(phi)
    pz = pavg*cos(theta)
    p.append(vector(px,py,pz)) 

# ---------------------------------------------------------
# HISTOGRAMA / DISTRIBUCIÓ D'ENERGIES.
# ---------------------------------------------------------
deltaE = (k*T) / 5 # Binning per l'histograma d'energies.
Emax = 6.0 * k * T # Límit màxim d'energies que es mostra a la gràfica.

def barE(E):
    return int(E/deltaE) # Funció que transforma una energía en un índex de l'histograma (similarment al que feia el codi original amb les velocitats).

n_bins = int(Emax/deltaE)
histo = []
for i in range(n_bins):
    histo.append(0.0)

Eavg_inicial = pavg**2 / (2*mass) # Introduim l'energia inicial (promig), a partir de pavg.

if barE(Eavg_inicial) < n_bins:
    histo[barE(Eavg_inicial)] = Natoms

gg = graph(width=win, height=0.4*win, xmax=Emax/1.6E-19, align='left', 
           xtitle='Energia (eV)', ytitle='N')
theory = gcurve(color=color.blue, width=2)


Etot_curve = gcurve(color=color.green)
t_sim = 0.0

# Predicció teòrica associada a la distribució de Maxwell Boltzmann per a Energies.
dE_step = Emax / 300.0
for i in range(301):
    E_val = i * dE_step
    if E_val == 0: E_val = 1e-30 # Evitamos división por 0
    # Fórmula teórica: f(E) = (2/sqrt(pi)) * (1/kT)^(3/2) * sqrt(E) * exp(-E/kT)
    f_E = (2.0/sqrt(pi)) * ((1.0/(k*T))**1.5) * sqrt(E_val) * exp(-E_val/(k*T))
    # Escalamos a número de partículas y tamaño de bin
    theory.plot(E_val/1.6E-19, Natoms * deltaE * f_E)

accum = []
for i in range(n_bins):
    accum.append([(deltaE*(i+0.5))/1.6E-19, 0])

# Gràfica de l'energia total (en eV) en funció del temps.
energy_graph = graph(width=win, height=0.4*win, align='left', 
                     xtitle='Passos de temps', ytitle='Energia total (eV)', ymin = 17.0, ymax = 24.0)
energy_curve = gcurve(color=color.green,width=2, graph=energy_graph)

# --- AÑADIDO: GRÁFICA DE FLUCTUACIONES RELATIVAS ---
fluct_graph = graph(width=win, height=0.4*win, align='left', 
                     xtitle='Passos de temps', ytitle='Fluctuacions rel. (%)')
fluct_curve = gcurve(color=color.orange, width=2, graph=fluct_graph)

# Calculamos la energía total esperada teórica
E_tot_esperada = Natoms * 1.5 * k * T
# ---------------------------------------------------

Edist = gvbars(color=color.red, delta=deltaE/1.6E-19, graph = gg)

def interchange(E1, E2): # Funció que actualitza l'histograma quan una partícula canvia d'energia.
    barE1 = barE(E1)
    barE2 = barE(E2)
    if barE1 == barE2: return
    if barE1 < len(histo):
        histo[barE1] -= 1
    if barE2 < len(histo):
        histo[barE2] += 1

def checkCollisions():
    hitlist = []
    r2 = 2*Ratom
    r2 *= r2
    for i in range(Natoms):
        ai = apos[i]
        for j in range(i) :
            aj = apos[j]
            dr = ai - aj
            if mag2(dr) < r2:
                hitlist.append([i,j])
    return hitlist

n_samples = 0 
t = 0 - dt #Inicialitzem la variable temporal.
nu = 5000 # Freqüència de col·lisió amb les partícules fictícies associades al bany tèrmic (és una mesura de l'acoblament amb el bany).
prob = nu * dt 

# Funció que genera una distribució gaussiana a partir d'una uniforme.
def box_muller():
    u1 = random.random()
    u2 = random.random()
    r = math.sqrt(-2 * math.log(u1))
    theta = 2 * math.pi * u2
    return r * math.cos(theta), r * math.sin(theta)

# BUCLE PRINCIPAL DE LA SIMULACIÓ.
while True:
    rate(300) 

    t += dt # Avancem en el temps.

    # 1) Actualitzem l'histograma.
    for i in range(len(accum)):
        accum[i][1] = (n_samples*accum[i][1] + histo[i])/(n_samples+1)
    if n_samples % 10 == 0:
        Edist.data = accum
    n_samples += 1

    E_total_julios = sum(mag2(p_i) / (2 * mass) for p_i in p)
    E_total_eV = E_total_julios / 1.6E-19
    energy_curve.plot(t/dt, E_total_eV)

    # --- AÑADIDO: PLOT DE LA FLUCTUACIÓN ---
    fluctuacio_percentatge = ((E_total_julios - E_tot_esperada) / E_tot_esperada) * 100
    fluct_curve.plot(t/dt, fluctuacio_percentatge)
    # ---------------------------------------

    # 2) Moviment de les partícules.
    for i in range(Natoms):
        Atoms[i].pos = apos[i] = apos[i] + (p[i]/mass)*dt

    
    
    # 3) Busquem quins àtoms han xocat.
    hitlist = checkCollisions()

    # 4) Si dos àtoms col·lisionen, acualitzem els seus moments.
    for ij in hitlist:
        i = ij[0]
        j = ij[1]
        ptot = p[i]+p[j]
        posi = apos[i]
        posj = apos[j]
        vi = p[i]/mass
        vj = p[j]/mass
        vrel = vj-vi
        a = vrel.mag2
        if a == 0: continue; 
        rrel = posi-posj
        if rrel.mag > Ratom: continue 

        # Guardem les energies inicials abans d'actualitzar moments.
        E_old_i = p[i].mag2 / (2*mass)
        E_old_j = p[j].mag2 / (2*mass)

        dx = dot(rrel, vrel.hat) 
        dy = cross(rrel, vrel.hat).mag 
        
        alpha = asin(dy/(2*Ratom)) 
        d = (2*Ratom)*cos(alpha)-dx 
        deltat = d/vrel.mag 
        
        posi = posi-vi*deltat 
        posj = posj-vj*deltat 
        mtot = 2*mass
        pcmi = p[i]-ptot*mass/mtot 
        pcmj = p[j]-ptot*mass/mtot 
        rrel = norm(rrel)
        pcmi = pcmi-2*pcmi.dot(rrel)*rrel 
        pcmj = pcmj-2*pcmj.dot(rrel)*rrel 
        p[i] = pcmi+ptot*mass/mtot 
        p[j] = pcmj+ptot*mass/mtot 
        apos[i] = posi+(p[i]/mass)*deltat 
        apos[j] = posj+(p[j]/mass)*deltat 

        # Calculem les noves energies i actualitzem l'histograma.
        E_new_i = p[i].mag2 / (2*mass)
        E_new_j = p[j].mag2 / (2*mass)
        interchange(E_old_i, E_new_i)
        interchange(E_old_j, E_new_j)

    for i in range(Natoms):
        loc = apos[i]
        if abs(loc.x) > L/2:
            if loc.x < 0: p[i].x = abs(p[i].x)
            else: p[i].x = -abs(p[i].x)
        if abs(loc.y) > L/2:
            if loc.y < 0: p[i].y = abs(p[i].y)
            else: p[i].y = -abs(p[i].y)
        if abs(loc.z) > L/2:
            if loc.z < 0: p[i].z = abs(p[i].z)
            else: p[i].z = -abs(p[i].z)

    # --------------------------------------
    # IMPLEMENTACIÓ DEL TERMOSTAT D'ANDERSEN
    # --------------------------------------
    for i in range(Natoms):
        if random.random() < prob: 
            # Guardem l'energía antiga per l'historial.
            E_old = p[i].mag2 / (2*mass)

            z1, z2 = box_muller()
            z3, _ = box_muller()
            
            sigma = math.sqrt(k*T/mass)
            
            vx = sigma * z1
            vy = sigma * z2
            vz = sigma * z3
            p[i] = vector(mass*vx, mass*vy, mass*vz)

            # Actualitzem l'histograma amb l'energía nova.
            E_new = p[i].mag2 / (2*mass)
            interchange(E_old, E_new)