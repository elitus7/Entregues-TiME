from vpython import *

win = 500

Natoms = 10  # Nombre d'àtoms ajustat a 500.

# Paràmetres de la simulació (tots els valors estan en SI).
L = 1 # Aresta de la caixa cúbica.
gray = color.gray(0.7) # Color de les arestes de la caixa.
mass = 4E-3/6E23 # Massa (en kg) d'un àtom d'He.
Ratom = 0.04 # Radi atòmic usat a la simulació. L'original és 0.03, el radi real és 0.31*10**-10 i el seleccionat per les simulacions és 0.04.
k = 1.4E-23 # Constant de Boltzmann.
T = 300 # Temperatura ambient (aproximadament).
dt = 1E-5

animation = canvas( width=win, height=win, align='left')
animation.range = L
animation.title = 'Hard Sphere Gas'
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
p = [] # Moments dels diferents àtoms.
apos = [] # Posicions dels diferents àtoms.
pavg = sqrt(2*mass*1.5*k*T) # Moment lineal associat a l'energia cinètica mitjana, p**2/(2mass) = (3/2)kT, que, recordem, queda fixada en seleccionar una T de treball (col·lectivitat NVE).
    
# Aquest bucle serveix per crear els diferents àtoms de la simulació, tot donant-lis una POSICIÓ inicial, una DIRECCIÓ en la que es mouran (en coordenades esfèriques) i el MOMENT inicial (que será el mateix per tots i es correspondrà amb pavg).
for i in range(Natoms):
    x = L*random()-L/2
    y = L*random()-L/2
    z = L*random()-L/2
    if i == 0:
        Atoms.append(sphere(pos=vector(x,y,z), radius=Ratom, color=color.cyan, make_trail=True, retain=100, trail_radius=0.3*Ratom)) #Al primer de tots els àtoms generats li posem una estela blava per poder seguir la seva posició (als altres no).
    else: Atoms.append(sphere(pos=vector(x,y,z), radius=Ratom, color=gray))
    apos.append(vec(x,y,z)) #Introduim les posicions donades a cada àtom (de forma aleatòria) a la llista de posicions [r_0, ..., r_N].
    theta = pi*random() #Donem un vector director aletaori a cada àtom (en funció dels angles theta i phi de les coordenades esfèriques).
    phi = 2*pi*random()
    px = pavg*sin(theta)*cos(phi) #A partir del vector director, donem les diferents components del vector moment lineal (tenint present que el seu módul és pavg).
    py = pavg*sin(theta)*sin(phi)
    pz = pavg*cos(theta)
    p.append(vector(px,py,pz)) #Incloem aquest moment a la llista de moments [p_0, ..., p_N].

# Aquesta part genera l'histograma de les velocitats de les partícules. Es representa també la predicció teórica associada a la distribució de Maxwell-Boltzmann.
deltav = 100 # Salt de velocitats a l'histograma.

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
for v in range(0,3001+dv,dv):  # Predicció teòrica associada a la distribució de Maxwell Boltzmann.
    theory.plot( v, (deltav/dv)*Natoms*4*pi*((mass/(2*pi*k*T))**1.5) *exp(-0.5*mass*(v**2)/(k*T))*(v**2)*dv )

accum = []
for i in range(int(3000/deltav)): accum.append([deltav*(i+.5),0])
vdist = gvbars(color=color.red, delta=deltav )

def interchange(v1, v2):  # Funció que actualitza l'histograma quan una partícula canvia de velocitat després d'una col·lisió.
    barx1 = barx(v1)
    barx2 = barx(v2)
    if barx1 == barx2:  return
    if barx1 >= len(histo) or barx2 >= len(histo): return
    histo[barx1] -= 1
    histo[barx2] += 1
    
def checkCollisions(): # Aquesta funció s'ocupa de determinar quines parelles d'àtoms col·lisionen
    hitlist = []
    r2 = 2*Ratom
    r2 *= r2
    for i in range(Natoms):
        ai = apos[i]
        for j in range(i) :
            aj = apos[j]
            dr = ai - aj
            if mag2(dr) < r2: hitlist.append([i,j]) # Notem que mag2(r) = x^2 + y^2 + z^2.
    return hitlist

nhisto = 0 # Nombre de cops que hem fet un histograma.
nu = 0.04 # Freqüència de col·lisió amb les partícules fictícies associades al bany tèrmic. És una mesura de l'acoblament entre el bany i el sistema.
prob = nu * dt # Probabilitat de que una partícula sigui seleccionada per interaccionar amb el bany tèrmic.

# BUCLE PRINCIPAL DE LA SIMULACIÓ.
while True:
    rate(300) # Fem que s'executi 300 cops per segon. Cada iteració farà el següent.

    # 1) Actualitzem l'histograma.
    for i in range(len(accum)): accum[i][1] = (nhisto*accum[i][1] + histo[i])/(nhisto+1)
    if nhisto % 10 == 0:
        vdist.data = accum
    nhisto += 1

    # 2) Fem que totes les partícules es moguin segons r_i = r_i(inicial) + v_i * dt = r_i(inicial) + p_i/m * dt.
    for i in range(Natoms): Atoms[i].pos = apos[i] = apos[i] + (p[i]/mass)*dt
    
    # 3) Busquem quins àtoms han xocat.
    hitlist = checkCollisions()

    # 4) Modifiquem el moment lineal de tots els àtoms que hagin col·lisionat, imposant que la col·lisió sigui ELÀSTICA. 
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
        if a == 0: continue;  # Velocitat relativa entre dues partícules és 0 --> tenen la mateixa velocitat. Si dues partícules amb la mateixa velocitat xoquen seguiran juntes sense que es produeixi cap "rebote".
        rrel = posi-posj
        if rrel.mag > Ratom: continue # Si la distància relativa és major que el radi atòmic, no hi ha col·lisió. Serveix per garantir que no hi hagi col·lisions falses.
    
        # L'angle theta es l'angle entre vrel i rrel. Si tenim això present, aleshores:
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
    
    for i in range(Natoms): # Aquesta part implementa els xocs amb les parets, fent que es produeixin reflexions totals.
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

    # --------------------------------------
    # IMPLEMENTACIÓ DEL TERMOSTAT D'ANDERSEN
    # --------------------------------------


    for i in range (Natoms):
        if random() < nu * dt:
            theta = pi*random()
            phi = 2*pi*random()

            p_mag = sqrt(2*mass*1.5*k*T)

            p[i] = vector(
                p_mag*sin(theta)*cos(phi),
                p_mag*sin(theta)*sin(phi),
                p_mag*cos(theta)
            )
    