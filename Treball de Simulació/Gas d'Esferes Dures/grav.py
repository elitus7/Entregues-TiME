from vpython import *

win = 500

Natoms = 500  # Nombre d'àtoms ajustat a 500.

# Paràmetres de la simulació (tots els valors estan en SI).
L = 1 # Aresta de la caixa cúbica.
gray = color.gray(0.7) # Color de les arestes de la caixa.
mass = 4E-3/6E23 # Massa (en kg) d'un àtom d'He.
Ratom = 0.04 # Radi atòmic usat a la simulació.
k = 1.4E-23 # Constant de Boltzmann.
T = 300 # Temperatura ambient (aproximadament).
g = 2000000 # Camp gravitatori (exagerat).
dt = 1E-5 # Pas de temps.


# Animació del sistema.
animation = canvas( width=win, height=win, align='left')
animation.range = L
animation.title = 'Hard Sphere Gas'
s = """  Distribució teòrica i simulada per les velocitats
  (en m/s) i posicions (en m) en l'eix Z amb camp gravitatori. 
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
    vz = pz/mass
    p.append(vector(px,py,pz))

# -----------------------------------
# HISTOGRAMA DE VELOCITATS EN L'EIX Z
# -----------------------------------
deltav = 200 # Amplada de les barres de l'histograma.

def barx(v):
    return int((v+3000)//deltav) 

nhisto_bars = int(6000/deltav)
histo = []
for i in range(nhisto_bars): histo.append(0.0)

for mom in p:
    idx = barx(mom.z/mass)
    if 0 <= idx < len(histo):
        histo[idx] += 1

# Creem l'entorn en què es mostrarà l'histograma en l'eix Z.
gg = graph( width=win, height=0.4*win, xmax=3000, align='left',
    xtitle='Velocitat eix Z (m/s)', ytitle='N', ymax=Natoms*deltav/1000, background=color.white)

theory = gcurve( color=color.blue, width=2 )
dv = 10
for vz in range(-3001+dv,3001+dv,dv):  
    theory.plot( vz, (deltav/dv)*Natoms*((mass/(2*pi*k*T))**0.5)*exp(-0.5*mass*(vz**2)/(k*T))*dv )

accum = []
for i in range(int(6000/deltav)): accum.append([-3000+deltav*(i+.5), 0])
vdist = gvbars(color=color.red, delta=deltav )

# ----------------------------------
# HISTOGRAMA DE POSICIONS EN L'EIX Z
# ----------------------------------
deltaz = 0.04 # Amplada de les barres de l'histograma (la caixa va de -0.5 a 0.5).
nhisto_z = int(1.0/deltaz)

def barz(z):
    return int((z + 0.5)/deltaz)

# Creem l'entorn en el que es generarà l'histograma de posicions en Z.
gg_pos = graph(width=win, height=0.4*win, xmax=0.5, xmin=-0.5, align='left',
    xtitle='Posició eix Z (m)', ytitle='N', ymax=50, background=color.white) 

accum_pos = []
for i in range(nhisto_z): accum_pos.append([-0.5 + deltaz*(i+.5), 0])
posdist = gvbars(color=color.orange, delta=deltaz)

# Funció que actualitza l'histograma quan una partícula canvia de velocitat. 
# Respecte al codi original això s'ha hagut de modificar lleugerament
# per evitar un descompte en el nombre de partícules de la simulació.
def interchange(v1, v2):  
    barx1 = barx(v1)
    barx2 = barx(v2)
    
    if barx1 == barx2: return
    
    if 0 <= barx1 < len(histo): 
        histo[barx1] -= 1
        
    if 0 <= barx2 < len(histo): 
        histo[barx2] += 1

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

nhisto = 0 
t = 0 - dt # Inicialitem variable temporal.

while True:
    rate(300)
    
    t += dt # Avancem en el temps un pas.

    # 1) Actualitzem histogrames.
    histo_pos = [0.0] * nhisto_z
    for loc in apos:
        histo_pos[barz(loc.z)] += 1

    for i in range(len(accum)): 
        accum[i][1] = (nhisto*accum[i][1] + histo[i])/(nhisto+1)
        
    for i in range(len(accum_pos)): 
        accum_pos[i][1] = (nhisto*accum_pos[i][1] + histo_pos[i])/(nhisto+1)

    if nhisto % 10 == 0:
        vdist.data = accum
        posdist.data = accum_pos 
    nhisto += 1

    # 2) Introduïm el moviment de les partícules i la gravetat.
    for i in range(Natoms): 
        # Guardem la component Z de la velocitat anterior per tal de posar-la a l'histograma.
        old_vz = p[i].z / mass  
        
        # Modifiquem la component Z del moment lineal de les partícules per introduir els efectes del camp gravitatori.
        p[i].z -= mass * g * dt 
        
        # Usant el nou moment, actualitzem les posicions.
        apos[i] = apos[i] + (p[i]/mass)*dt
        Atoms[i].pos = apos[i]
        
        # Actualitzem l'histograma.
        interchange(old_vz, p[i].z / mass)
    
    # 3) Busquem quines partícules col·lisionen.
    hitlist = checkCollisions()

    # 4) Si dos àtoms col·lisionen, actualitzem les seves posicions i moments.
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
        
        interchange(vi.z, p[i].z/mass)
        interchange(vj.z, p[j].z/mass)
    
    # 5) Col·lisions amb les parets (part modificada per garantir millor estabilitat numèrica en introduir el camp gravitatori).
    for i in range(Natoms):
        loc = apos[i]
        
        # Rebot en X.
        if abs(loc.x) > L/2:
            if loc.x < 0: 
                p[i].x = abs(p[i].x)
                apos[i].x = -L - loc.x 
            else: 
                p[i].x = -abs(p[i].x)
                apos[i].x = L - loc.x  
            Atoms[i].pos.x = apos[i].x
        
        # Rebot en Y.
        if abs(loc.y) > L/2:
            if loc.y < 0: 
                p[i].y = abs(p[i].y)
                apos[i].y = -L - loc.y 
            else: 
                p[i].y = -abs(p[i].y)
                apos[i].y = L - loc.y  
            Atoms[i].pos.y = apos[i].y
        
        # Rebot en Z.
        if abs(loc.z) > L/2:
            old_vz = p[i].z / mass
            if loc.z < 0: 
                p[i].z = abs(p[i].z)
                apos[i].z = -L - loc.z 
            else: 
                p[i].z = -abs(p[i].z)
                apos[i].z = L - loc.z  
            Atoms[i].pos.z = apos[i].z
            
            interchange(old_vz, p[i].z / mass)