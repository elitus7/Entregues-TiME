from vpython import *

win = 500
Natoms = 500 # Nombre d'àtoms ajustat a 500.

# Paràmetres de la simulació (tots els valors estan en SI).
L = 1 # Aresta de la caixa cúbica.
gray = color.gray(0.7) # Color de les arestes de la caixa.
mass = 4E-3/6E23 # Massa (en kg) d'un àtom d'He.
Ratom = 0.04 # Radi atòmic usat a la simulació.
k = 1.4E-23 # Constant de Boltzmann.
T = 300 # Temperatura ambient (aproximadament).
dt = 1E-5

# --- PARÁMETRO NUEVO: TERMOSTATO DE ANDERSEN ---
nu = 10000 # Frecuencia de colisión con el baño térmico (ajustable).
# -----------------------------------------------

animation = canvas( width=win, height=win, align='left')
animation.range = L
animation.title = 'Hard Sphere Gas with Andersen Thermostat'
s = """  Theoretical and averaged speed distributions (meters/sec).
  Initially all atoms have the same speed, but collisions
  change the speeds of the colliding atoms. One of the atoms is
  marked and leaves a trail so you can follow its path.
  
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
pavg = sqrt(2*mass*1.5*k*T) # average kinetic energy p**2/(2mass) = (3/2)kT
    
for i in range(Natoms):
    x = L*random()-L/2
    y = L*random()-L/2
    z = L*random()-L/2
    if i == 0:
        Atoms.append(sphere(pos=vector(x,y,z), radius=Ratom, color=color.cyan, make_trail=True, retain=100, trail_radius=0.3*Ratom))
    else: Atoms.append(sphere(pos=vector(x,y,z), radius=Ratom, color=gray))
    apos.append(vec(x,y,z))
    theta = pi*random()
    phi = 2*pi*random()
    px = pavg*sin(theta)*cos(phi)
    py = pavg*sin(theta)*sin(phi)
    pz = pavg*cos(theta)
    p.append(vector(px,py,pz))

deltav = 100 # binning for v histogram

def barx(v):
    return int(v/deltav) # index into bars array

nhisto = int(4500/deltav)
histo = []
for i in range(nhisto): histo.append(0.0)
histo[barx(pavg/mass)] = Natoms

gg = graph( width=win, height=0.4*win, xmax=3000, align='left',
    xtitle='Velocitat (m/s)', ytitle='N', ymax=Natoms*deltav/1000)

theory = gcurve( color=color.blue, width=2 )
dv = 10
for v in range(0,3001+dv,dv):  # theoretical prediction
    theory.plot( v, (deltav/dv)*Natoms*4*pi*((mass/(2*pi*k*T))**1.5) *exp(-0.5*mass*(v**2)/(k*T))*(v**2)*dv )

accum = []
for i in range(int(3000/deltav)): accum.append([deltav*(i+.5),0])
vdist = gvbars(color=color.red, delta=deltav )

def interchange(v1, v2):  # remove from v1 bar, add to v2 bar
    barx1 = barx(v1)
    barx2 = barx(v2)
    if barx1 == barx2:  return
    if barx1 >= len(histo) or barx2 >= len(histo): return
    histo[barx1] -= 1
    histo[barx2] += 1
    
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

# --- FUNCIÓN NUEVA: GENERADOR MAXWELL-BOLTZMANN ---
def get_MB_velocity(m, T_desired):
    # Desviación estándar de la distribución de velocidades
    sigma = sqrt(k * T_desired / m)
    
    # Transformada de Box-Muller para V_x
    u1 = random(); 
    while u1 == 0: u1 = random() # Evitar log(0)
    u2 = random()
    vx = sigma * sqrt(-2*log(u1)) * cos(2*pi*u2)
    
    # Transformada de Box-Muller para V_y
    u1 = random(); 
    while u1 == 0: u1 = random()
    u2 = random()
    vy = sigma * sqrt(-2*log(u1)) * cos(2*pi*u2)
    
    # Transformada de Box-Muller para V_z
    u1 = random(); 
    while u1 == 0: u1 = random()
    u2 = random()
    vz = sigma * sqrt(-2*log(u1)) * cos(2*pi*u2)
    
    return vector(vx, vy, vz)
# --------------------------------------------------

nhisto = 0 # number of histogram snapshots to average

while True:
    rate(300)
    # Accumulate and average histogram snapshots
    for i in range(len(accum)): accum[i][1] = (nhisto*accum[i][1] + histo[i])/(nhisto+1)
    if nhisto % 10 == 0:
        vdist.data = accum
    nhisto += 1

    # Update all positions (Paso 1 del termostato de Andersen)
    for i in range(Natoms): Atoms[i].pos = apos[i] = apos[i] + (p[i]/mass)*dt
    
    # Check for collisions
    hitlist = checkCollisions()

    # If any collisions took place, update momenta of the two atoms
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
        interchange(vi.mag, p[i].mag/mass)
        interchange(vj.mag, p[j].mag/mass)
    
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

    # --- IMPLEMENTACIÓN DEL TERMOSTATO DE ANDERSEN ---
    # (Pasos 2 y 3 de tu imagen)
    prob_colision = nu * dt 
    
    for i in range(Natoms):
        if random() < prob_colision:
            v_old_mag = p[i].mag / mass
            
            # Obtener nueva velocidad y actualizar momento
            v_new = get_MB_velocity(mass, T)
            p[i] = v_new * mass
            v_new_mag = p[i].mag / mass
            
            # Actualizar el histograma visual para reflejar el cambio
            interchange(v_old_mag, v_new_mag)
    # -------------------------------------------------
    # Añade esto al final del bucle while True:
    if nhisto % 100 == 0: # Imprimir cada 100 ciclos para no saturar la consola
        K_total = 0
        for i in range(Natoms):
            K_total += p[i].mag2 / (2 * mass)
        
        T_inst = (2/3) * K_total / (Natoms * k)
        print(f"Temperatura instantánea: {T_inst:.2f} K")