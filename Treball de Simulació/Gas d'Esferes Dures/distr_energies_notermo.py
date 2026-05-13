from vpython import *
import numpy as np
import math
import random

win = 500

Natoms = 500  # Nombre d'àtoms ajustat a 500.

# Paràmetres de la simulació (tots els valors estan en SI).
L = 1 # Aresta de la caixa cúbica.
gray = color.gray(0.7) # Color de les arestes de la caixa.
mass = 4E-3/6E23 # Massa (en kg) d'un àtom d'He.
Ratom = 0.04 # Radi atòmic usat a la simulació.
k = 1.4E-23 # Constant de Boltzmann.
T = 300 # Temperatura ambient (aproximadament).
dt = 1E-5 # Pas de temps.

# Animació del sistema.
animation = canvas( width=win, height=win, align='left')
animation.range = L
animation.title = 'Hard Sphere Gas'
s = """  Distribució teòrica i simulada per les energies
  (en eV), energia en funció del temps i fluctiacions relatives 
  en funció del temps. Inicialment tots els àtoms tenen la mateixa
  velocitat, però les col·lisions fan que canviïn. Un dels àtoms
  està marcat per tal de poder seguir la seva trajectòria.

  Col·lectivitat NVE.
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

Atoms = [] # Llista amb tots els àtoms.
p = [] # Llista de tots els moments dels àtoms.
apos = [] # Llista de totes les velocitats dels àtoms.
pavg = sqrt(2*mass*1.5*k*T) # Energia cinètica promig: p**2/(2mass) = (3/2)kT.

for i in range(Natoms):
    x = L*random.random()-L/2
    y = L*random.random()-L/2
    z = L*random.random()-L/2
    if i == 0:
        Atoms.append(sphere(pos=vector(x,y,z), radius=Ratom, color=color.cyan, make_trail=True, retain=100, trail_radius=0.3*Ratom))
    else: Atoms.append(sphere(pos=vector(x,y,z), radius=Ratom, color=gray))
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

# Predicció teòrica associada a la distribució de Maxwell-Boltzmann per energies.
dE_step = Emax / 300.0
for i in range(301):
    E_val = i * dE_step
    f_E = (2.0/sqrt(pi)) * ((1.0/(k*T))**1.5) * sqrt(E_val) * exp(-E_val/(k*T))
    theory.plot(E_val/1.6E-19, Natoms * deltaE * f_E)

accum = []
for i in range(n_bins):
    accum.append([(deltaE*(i+0.5))/1.6E-19, 0])

# Gràfica de l'energia total (en eV) en funció del temps.
energy_graph = graph(width=win, height=0.4*win, align='left', 
                     xtitle='Passos de temps', ytitle='Energia total (eV)', ymin = 19.5, ymax = 19.8)
energy_curve = gcurve(color=color.green, width=2, graph=energy_graph)

Edist = gvbars(color=color.red, delta=deltaE/1.6E-19, graph = gg)

# Funció que actualitza l'histograma quan una partícula canvia d'energia.
def interchange(E1, E2): 
    barE1 = barE(E1)
    barE2 = barE(E2)
    if barE1 == barE2: return
    if barE1 < len(histo):
        histo[barE1] -= 1
    if barE2 < len(histo):
        histo[barE2] += 1

# Funció que comprova si dues partícules han col·lisionat.
def checkCollisions(): 
    hitlist = []
    r2 = 2*Ratom
    r2 *= r2
    for i in range(Natoms):
        ai = apos[i]
        for j in range(i) :
            aj = apos[j]
            dr = ai - aj
            if mag2(dr) < r2: hitlist.append([i,j])
    return hitlist

n_samples = 0 
t = 0 - dt # Inicialitzem la variable temporal.

# -------------------------------
# BUCLE PRINCIPAL DE LA SIMULACIÓ
# -------------------------------
while True:
    rate(300)
    
    t += dt # Avancem en el temps un pas.
    
    # 1) Actualitzem l'histograma i calculem l'energia total (en eV). Fem plot al gràfic E(t)
    for i in range(len(accum)):
        accum[i][1] = (n_samples*accum[i][1] + histo[i])/(n_samples+1)
    if n_samples % 10 == 0:
        Edist.data = accum
    n_samples += 1
    
    E_total_jules = sum(mag2(p_i) / (2 * mass) for p_i in p)
    E_total_eV = E_total_jules / 1.6E-19
    energy_curve.plot(t/dt, E_total_eV)

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
        if a == 0: continue;  # exactly same velocities
        rrel = posi-posj
        if rrel.mag > Ratom: continue # one atom went all the way through another

        # Guardem les energies inicials abans d'actualitzar moments.
        E_old_i = p[i].mag2 / (2*mass)
        E_old_j = p[j].mag2 / (2*mass)

        # Calculem les noves posicions i moments després de la col·lisió.
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
            if loc.x < 0: p[i].x =  abs(p[i].x)
            else: p[i].x =  -abs(p[i].x)
        
        if abs(loc.y) > L/2:
            if loc.y < 0: p[i].y = abs(p[i].y)
            else: p[i].y =  -abs(p[i].y)
        
        if abs(loc.z) > L/2:
            if loc.z < 0: p[i].z =  abs(p[i].z)
            else: p[i].z =  -abs(p[i].z)