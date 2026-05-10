from vpython import *

win = 500

Natoms = 500  # Nombre d'àtoms ajustat a 500.

# Paràmetres de la simulació (tots els valors estan en SI).
L = 1 # Aresta de la caixa cúbica.
gray = color.gray(0.7) # Color de les arestes de la caixa.
mass = 4E-3/6E23 # Massa (en kg) d'un àtom d'He.
Ratom = 0.04 # Radi atòmic usat a la simulació. L'original és 0.03, el radi real és 0.31*10**-10 i el seleccionat per les simulacions és 0.04.
k = 1.4E-23 # Constant de Boltzmann.
T = 300 # Temperatura ambient (aproximadament).
dt = 1E-5 # Pas de temps.

animation = canvas( width=win, height=win, align='left')
animation.range = L
animation.title = 'Hard Sphere Gas'
s = """  Distribució teòrica i simulada per les velocitats
  (en m/s). Inicialment tots els àtoms tenen la mateixa
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

deltav = 100 # Binning per l'histograma d'energies.

def barx(v):
    return int(v/deltav) # Funció que transforma una velocitat en un índex de l'histograma.

nhisto = int(4500/deltav)
histo = []
for i in range(nhisto): histo.append(0.0)
histo[barx(pavg/mass)] = Natoms

gg = graph( width=win, height=0.4*win, xmax=3000, align='left',
    xtitle='Velocitat (m/s)', ytitle='N', ymax=Natoms*deltav/1000)

theory = gcurve( color=color.blue, width=2 )
dv = 10
for v in range(0,3001+dv,dv):  # Predicció teòrica.
    theory.plot( v, (deltav/dv)*Natoms*4*pi*((mass/(2*pi*k*T))**1.5) *exp(-0.5*mass*(v**2)/(k*T))*(v**2)*dv )

accum = []
for i in range(int(3000/deltav)): accum.append([deltav*(i+.5),0])
vdist = gvbars(color=color.red, delta=deltav )

def interchange(v1, v2):  # # Funció que actualitza l'histograma quan una partícula canvia d'energia.
    barx1 = barx(v1)
    barx2 = barx(v2)
    if barx1 == barx2:  return
    if barx1 >= len(histo) or barx2 >= len(histo): return
    histo[barx1] -= 1
    histo[barx2] += 1
    
def checkCollisions(): # Funció que comprova si dues partícules han col·lisionat.
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

nhisto = 0 

# -------------------------------
# BUCLE PRINCIPAL DE LA SIMULACIÓ
# -------------------------------
while True:
    rate(300)
    # # 1) Actualitzem l'histograma i calculem l'energia total (en eV). Fem els plots corresponents als gràfics.
    for i in range(len(accum)): accum[i][1] = (nhisto*accum[i][1] + histo[i])/(nhisto+1)
    if nhisto % 10 == 0:
        vdist.data = accum
    nhisto += 1

    # 2) Moviment de les partícules.
    for i in range(Natoms): Atoms[i].pos = apos[i] = apos[i] + (p[i]/mass)*dt
    
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
    
        # theta is the angle between vrel and rrel:
        dx = dot(rrel, vrel.hat)       # rrel.mag*cos(theta)
        dy = cross(rrel, vrel.hat).mag # rrel.mag*sin(theta)
        # alpha is the angle of the triangle composed of rrel, path of atom j, and a line
        #   from the center of atom i to the center of atom j where atome j hits atom i:
        alpha = asin(dy/(2*Ratom)) 
        d = (2*Ratom)*cos(alpha)-dx # distance traveled into the atom from first contact
        deltat = d/vrel.mag         # time spent moving from first contact to position inside atom
        
        posi = posi-vi*deltat # back up to contact configuration
        posj = posj-vj*deltat
        mtot = 2*mass
        pcmi = p[i]-ptot*mass/mtot # transform momenta to cm frame
        pcmj = p[j]-ptot*mass/mtot
        rrel = norm(rrel)
        pcmi = pcmi-2*pcmi.dot(rrel)*rrel # bounce in cm frame
        pcmj = pcmj-2*pcmj.dot(rrel)*rrel
        p[i] = pcmi+ptot*mass/mtot # transform momenta back to lab frame
        p[j] = pcmj+ptot*mass/mtot
        apos[i] = posi+(p[i]/mass)*deltat # move forward deltat in time
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

    


