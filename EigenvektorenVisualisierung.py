#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Fortgeschrittenen Praktikum TU Berlin
# P60 Ramanspektroskopie an Halbleitern
# WS19/20 Gaebel & Krueger
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):

    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)

def sphere(position,radius,col):    
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = radius*np.cos(u)*np.sin(v)
    y = radius*np.sin(u)*np.sin(v)
    z = radius*np.cos(v)
    ax.plot_wireframe(x+position[0], y+position[1], z+position[2], color=col)

#-- Silizium --#
fcc1 = np.array([
        [0,0,0],
        [1,0,0],
        [0,1,0],
        [0,0,1],
        
        [1,1,0],
        [0,1,1],
        [1,0,1],
        [1,1,1],
        
        [0.5,0.5,0],
        [0.5,0,0.5],
        [0,0.5,0.5],
        [1,0.5,0.5],
        [0.5,1,0.5],
        [0.5,0.5,1]
        ])

fcc2 = np.array([
        [0.25,0.25,0.25],
        [0.75,0.75,0.25],
        [0.75,0.25,0.75],
        [0.25,0.75,0.75],
        ])

border = np.array([
        [[0,0,0],[1,0,0]],
        [[1,0,0],[1,1,0]],
        [[1,1,0],[0,1,0]],
        [[0,1,0],[0,0,0]],
        
        [[0,0,0],[0,0,1]],
        [[1,0,0],[1,0,1]],
        [[1,1,0],[1,1,1]],
        [[0,1,0],[0,1,1]],
        
        [[0,0,1],[1,0,1]],
        [[1,0,1],[1,1,1]],
        [[1,1,1],[0,1,1]],
        [[0,1,1],[0,0,1]],
        ])

connections = np.array([
        [[0.25,0.25,0.25],[0,0,0]],   
        [[0.25,0.25,0.25],[0.5,0.5,0]],
        [[0.25,0.25,0.25],[0.5,0,0.5]],
        [[0.25,0.25,0.25],[0,0.5,0.5]],
        
        [[0.75,0.75,0.25],[0.5,0.5,0]],
        [[0.75,0.75,0.25],[1,0.5,0.5]],
        [[0.75,0.75,0.25],[0.5,1,0.5]],
        [[0.75,0.75,0.25],[1,1,0]],
        
        [[0.75,0.25,0.75],[0.5,0,0.5]],
        [[0.75,0.25,0.75],[1,0.5,0.5]],
        [[0.75,0.25,0.75],[0.5,0.5,1]],
        [[0.75,0.25,0.75],[1,0,1]],
               
        [[0.25,0.75,0.75],[0,0.5,0.5]],
        [[0.25,0.75,0.75],[0.5,1,0.5]],
        [[0.25,0.75,0.75],[0.5,0.5,1]],
        [[0.25,0.75,0.75],[0,1,1]],
        ])
    
# Eigenvektoren der Phononen Aufschlüsseln
    
# Ausgabe:
'''         freq (    1) =       0.117738 [THz] =       3.927312 [cm-1]
 ( -0.520084  0.000000  0.479074 -0.000000 -0.000000  0.000000 ) 
 ( -0.520084  0.000000  0.479074 -0.000000 -0.000000  0.000000 ) 
     freq (    2) =       0.117738 [THz] =       3.927312 [cm-1]
 ( -0.000000  0.000000  0.000000  0.000000  0.707107  0.000000 ) 
 ( -0.000000  0.000000 -0.000000  0.000000  0.707107 -0.000000 ) 
     freq (    3) =       0.117738 [THz] =       3.927312 [cm-1]
 ( -0.479074  0.000000 -0.520084  0.000000 -0.000000  0.000000 ) 
 ( -0.479074  0.000000 -0.520084  0.000000 -0.000000  0.000000 ) 
     freq (    4) =      14.724612 [THz] =     491.160186 [cm-1]
 ( -0.302135  0.000000 -0.049077  0.000000 -0.637421  0.000000 ) 
 (  0.302135  0.000000  0.049077  0.000000  0.637421  0.000000 ) 
     freq (    5) =      14.724612 [THz] =     491.160186 [cm-1]
 ( -0.593744  0.000000 -0.239850  0.000000  0.299899 -0.000000 ) 
 (  0.593744 -0.000000  0.239850 -0.000000 -0.299899  0.000000 )
      freq (    6) =      14.724612 [THz] =     491.160186 [cm-1]
 (  0.237027  0.000000 -0.663373  0.000000 -0.061275  0.000000 ) 
 ( -0.237027  0.000000  0.663373  0.000000  0.061275  0.000000 )'''

# f in THz, f in 1/cm, v1, v2
eigenVectors=np.array([
        [0.117738,3.927312,[-0.520084,0.479074,0.0],[-0.520084,0.479074,0.0]],
        [0.117738,3.927312,[0.0,0.0,0.707107],[0.0,0.0,0.707107]],
        [0.117738,3.927312,[-0.479074,-0.520084,0.0],[-0.479074,-0.520084,0.0]],
        [14.724612,491.160186,[-0.302135,-0.049077,-0.637421],[0.302135,0.049077,0.637421]],
        [14.724612,491.160186,[-0.593744,-0.239850,0.299899],[0.593744,0.239850,-0.299899]],
        [14.724612,491.160186,[0.237027,-0.663373,-0.061275],[-0.237027,0.663373,0.061275]],
        ])

for numberEV in range(6):    
    fig = plt.figure()
    ax = fig.gca(projection='3d')    
    ax._axis3don = False
    
    # Plotten des Kristalls
    radius = 0.05
    for i in fcc1:
        sphere(i,radius,'g')
    for i in fcc2:
        sphere(i,radius,'b')
    
    for i in border:
        a = Arrow3D([i[0,0], i[1,0]], [i[0,1], i[1,1]], [i[0,2], i[1,2]], mutation_scale=20,
                    lw=1, arrowstyle="-", color="k")
        ax.add_artist(a)
    for i in connections:
        a = Arrow3D([i[0,0], i[1,0]], [i[0,1], i[1,1]], [i[0,2], i[1,2]], mutation_scale=20,
                    lw=3, arrowstyle="-", color="r")
        ax.add_artist(a)
    
    # Plotten der Eigenvektoren    
    normlength = 0.3 
    # Gitter-Punkte
    vec = eigenVectors[numberEV,2]
    # Normierung der EV
    normVec = (normlength/np.sqrt(vec[0]**2+vec[1]**2+vec[2]**2))*np.array(vec)
    for pos in fcc1:
        a = Arrow3D([pos[0], pos[0]+normVec[0]], [pos[1], pos[1]+normVec[1]], [pos[2], pos[2]+normVec[2]], mutation_scale=20,lw=2, arrowstyle="->", color="g")
        ax.add_artist(a)
    # Neben-Punkte
    vec = eigenVectors[numberEV,3]
    normVec = (normlength/np.sqrt(vec[0]**2+vec[1]**2+vec[2]**2))*np.array(vec)
    for pos in fcc2:
        a = Arrow3D([pos[0], pos[0]+normVec[0]], [pos[1], pos[1]+normVec[1]], [pos[2], pos[2]+normVec[2]], mutation_scale=20,lw=2, arrowstyle="->", color="b")
        ax.add_artist(a)
        
        plt.savefig('SiPhononenEV'+str(numberEV)+'.pdf')
        plt.show()

#-- Graphen --#

# Gittervektoren
e1 = np.array([1,0,0])
e2 = np.array([-0.5,np.sqrt(3)/2,0])

# Gitter Erzeugen
xdim = 3
ydim = 3
# Gitter1 generieren
hex1 = np.zeros([xdim,ydim,3])
for x in range(xdim):
    for y in range(ydim):
        hex1[x,y,:] = x*e1+y*e2
# Gitter2 gnerieren
hex2 = np.zeros([xdim,ydim,3])
for x in range(xdim):
    for y in range(ydim):
        hex2[x,y,:] = x*e1+y*e2+e1*1/3+e2*2/3
    
# Eigenvektoren der Phononen Aufschlüsseln
    
# Ausgabe:
"""     freq (    1) =       2.174422 [THz] =      72.530903 [cm-1]
 (  0.639581  0.000000  0.301556 -0.000000  0.000000 -0.000000 ) 
 (  0.639581  0.000000  0.301556 -0.000000  0.000000 -0.000000 ) 
     freq (    2) =       2.174422 [THz] =      72.530903 [cm-1]
 (  0.301556 -0.000000 -0.639581  0.000000 -0.000000  0.000000 ) 
 (  0.301556 -0.000000 -0.639581  0.000000 -0.000000  0.000000 ) 
     freq (    3) =       2.220315 [THz] =      74.061745 [cm-1]
 ( -0.000000  0.000000 -0.000000  0.000000  0.707107  0.000000 ) 
 ( -0.000000  0.000000 -0.000000  0.000000  0.707107  0.000000 ) 
     freq (    4) =      26.506295 [THz] =     884.154835 [cm-1]
 ( -0.000000  0.000000 -0.000000  0.000000  0.707107  0.000000 ) 
 (  0.000000  0.000000 -0.000000  0.000000 -0.707107  0.000000 ) 
     freq (    5) =      47.408535 [THz] =    1581.378509 [cm-1]
 (  0.002109  0.000000 -0.707104  0.000000 -0.000000  0.000000 ) 
 ( -0.002109  0.000000  0.707104 -0.000000  0.000000  0.000000 ) 
     freq (    6) =      47.408535 [THz] =    1581.378509 [cm-1]
 ( -0.707104  0.000000 -0.002109  0.000000 -0.000000  0.000000 ) 
 (  0.707104 -0.000000  0.002109  0.000000  0.000000 -0.000000 )"""

# f in THz, f in 1/cm, v1, v2
eigenVectors=np.array([
        [2.174422,72.530903,[0.639581,0.301556,0.000000],[0.639581,0.301556,0.000000]],
        [2.174422,72.530903,[0.301556,-0.639581,-0.000000],[0.301556,-0.639581,-0.000000]],
        [2.220315,74.061745,[-0.000000,-0.000000,0.707107],[-0.000000,-0.000000,0.707107]],
        [26.506295,884.154835,[-0.000000,-0.000000,0.707107],[0.000000,-0.000000,-0.707107]],
        [47.408535,1581.378509,[0.002109,-0.707104,-0.000000],[-0.002109,0.707104,0.000000]],
        [47.408535,1581.378509,[-0.707104,-0.002109,-0.000000],[0.707104,0.002109,0.000000]],
        ])
  
for numberEV in range(6):
    fig = plt.figure(constrained_layout=True)
    ax = fig.gca(projection='3d')   
    ax.set_zlim(-1, 1)
    ax._axis3don = False    
    
    # Plotten des Kristalls
    # Atome
    radius = 0.05
    for i in hex1:
        for j in i:
            sphere(j,radius,'g')
    for i in hex2:
        for j in i:
            sphere(j,radius,'b')      
    # Verbindungen
    for i in hex1:
        for j in i:
            a = Arrow3D([j[0],j[0]],[j[1],j[1]+1/np.sqrt(3)],[j[2],j[2]], mutation_scale=20,
                        lw=1, arrowstyle="-", color="k")
            ax.add_artist(a)
            v = e2-[0,1/np.sqrt(3),0]
            a = Arrow3D([j[0],j[0]+v[0]],[j[1],j[1]-v[1]],[j[2],j[2]+v[2]], mutation_scale=20,
                lw=1, arrowstyle="-", color="k")
            ax.add_artist(a)
            a = Arrow3D([j[0],j[0]-v[0]],[j[1],j[1]-v[1]],[j[2],j[2]+v[2]], mutation_scale=20,
                lw=1, arrowstyle="-", color="k")
            ax.add_artist(a)        
    for i in hex2:
        for j in i:
            a = Arrow3D([j[0],j[0]],[j[1],j[1]-1/np.sqrt(3)],[j[2],j[2]], mutation_scale=20,
                        lw=1, arrowstyle="-", color="k")
            ax.add_artist(a)
            v = e2-[0,1/np.sqrt(3),0]
            a = Arrow3D([j[0],j[0]-v[0]],[j[1],j[1]+v[1]],[j[2],j[2]+v[2]], mutation_scale=20,
                lw=1, arrowstyle="-", color="k")
            ax.add_artist(a)
            a = Arrow3D([j[0],j[0]+v[0]],[j[1],j[1]+v[1]],[j[2],j[2]+v[2]], mutation_scale=20,
                lw=1, arrowstyle="-", color="k")
            ax.add_artist(a)
                
    # Plotten der Eigenvektoren    
    normlength = 0.4
    # Gitter-Punkte
    vec = eigenVectors[numberEV,2]
    # Normierung der EV
    normVec = (normlength/np.sqrt(vec[0]**2+vec[1]**2+vec[2]**2))*np.array(vec)
    for x in range(xdim):
        for y in range(ydim):
            pos = hex1[x,y,:]
            a = Arrow3D([pos[0], pos[0]+normVec[0]], [pos[1], pos[1]+normVec[1]], [pos[2], pos[2]+normVec[2]], mutation_scale=20,lw=2, arrowstyle="->", color="g")
            ax.add_artist(a)
    # Neben-Punkte
    vec = eigenVectors[numberEV,3]
    normVec = (normlength/np.sqrt(vec[0]**2+vec[1]**2+vec[2]**2))*np.array(vec)
    for x in range(xdim):
        for y in range(ydim):
            pos = hex2[x,y,:]
            a = Arrow3D([pos[0], pos[0]+normVec[0]], [pos[1], pos[1]+normVec[1]], [pos[2], pos[2]+normVec[2]], mutation_scale=20,lw=2, arrowstyle="->", color="b")
            ax.add_artist(a)
    
    plt.savefig('CPhononenEV'+str(numberEV)+'.pdf')
    plt.show()