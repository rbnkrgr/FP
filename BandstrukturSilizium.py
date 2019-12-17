#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Fortgeschrittenen Praktikum TU Berlin
# P60 Ramanspektroskopie an Halbleitern
# WS19/20 Gaebel & Krueger
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc

## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

a0 = 0.529177210903 # Borscher radius in Angström

def floatFromString(string):
    '''Die in einem String mit anderen Buchstaben enthaltene Zahl wird
    herausgesucht und in eine float konvertiert.'''
    beginningList = ['0','1','2','3','4','5','6','7','8','9','-']
    continueList = ['0','1','2','3','4','5','6','7','8','9','.','-']
    for i in range(len(string)):
        if string[i] in beginningList:
            index1 = i
            break
        i = i+1
    for i in range(index1,len(string)):
        if not (string[i] in continueList):
            index2 = i
            break
        i = i+1
    return float(string[index1:index2])  

# Pfade werden festgelegt
pathQE = 'pw.x'

# Inputdatei wird generiert
stringIN_scf = '''      
        &control
        calculation = 'scf',
        ! Berechnungsart: scf-Berechnung wird durchgeführt
     /
     &system
        ibrav = 2,
        ! Gitter-Struktur: fcc (diamant ist fcc mit 2 Basisatomen)
        celldm(1) = 10.4,
        ! Gitterkonstante a in Bohr wird fetgelegt
        nat = 2,
        ! Anzahl an Atomen pro Elementarzelle
        ntyp= 1,
        ! Anzahl an Atomsorten
        ecutwfc = 40,
        ! Energie, bis zu welcher die Funktionen entwickelt werden in Ry (CutOff)
	ecutrho = 320,
	! maximale kinetische Energie der Elektronen (hier 8 x ecutwfc)
     /
     &electrons
        mixing_beta = 0.7
        ! Gewichtung zwischen aufeinanderfolgenden Dichteverteilungen
     /
    ATOMIC_SPECIES
     Si  28.086 Si.rel-pbe-rrkj.UPF
     ! Element, Masse in a.m.u., Pseudopotential
    ATOMIC_POSITIONS alat
     Si 0.00 0.00 0.00
     Si 0.25 0.25 0.25   
     ! alat: Positionen werden in relativen Einheiten der Gitterkonstante angegeben
    K_POINTS crystal_b
     5
	L 20
	gG 30
	X 10
	U 30
	gG 20
 '''

# Eingabedatei wird erstellt   
nameIN_scf = 'Input_Si_scf.in'
file = open(nameIN_scf,'w')
file.write(stringIN_scf)
file.close()  

# Name der Ausgabedatei 
nameOUT_scf = 'Ouput_Si_scf.out'    

# Simulation wird durchgeführt
os.system('pw.x <'+nameIN_scf+' >'+nameOUT_scf)


stringIN_bands = '''      
            &control
        calculation = 'bands',
        ! Berechnungsart: Bandstruktur wir berechnet
     /
     &system
        ibrav = 2,
        ! Gitter-Struktur: fcc (diamant ist fcc mit 2 Basisatomen)
        celldm(1) = 10.4,
        ! Gitterkonstante a in Bohr wird fetgelegt
        nat = 2,
        ! Anzahl an Atomen pro Elementarzelle
        ntyp= 1,
        ! Anzahl an Atomsorten
        ecutwfc = 40,
        ! Energie, bis zu welcher die Funktionen entwickelt werden in Ry (CutOff)
	ecutrho = 320,
	! maximale kinetische Energie der Elektronen (hier 8 x ecutwfc)
	nbnd = 8,
	! Anzahl der berechneten Bänder
     /
     &electrons
        mixing_beta = 0.7
        ! Gewichtung zwischen aufeinanderfolgenden Dichteverteilungen
	!conv_thr = 1e-8
     /
    ATOMIC_SPECIES
     Si  28.086 Si.rel-pbe-rrkj.UPF
     ! Element, Masse in a.m.u., Pseudopotential
    ATOMIC_POSITIONS alat
     Si 0.00 0.00 0.00
     Si 0.25 0.25 0.25   
     ! alat: Positionen werden in relativen Einheiten der Gitterkonstante angegeben
    K_POINTS crystal_b
     5
	L 20
	gG 30
	X 10
	U 30
	gG 20
 '''

# Eingabedatei wird erstellt   
nameIN_bands = 'Input_Si_bands.in'
file = open(nameIN_bands,'w')
file.write(stringIN_bands)
file.close()  

# Name der Ausgabedatei 
nameOUT_bands = 'Ouput_Si_bands.out'    

# Simulation wird durchgeführt
os.system('pw.x <'+nameIN_bands+' >'+nameOUT_bands)


stringIN_bands_mod = '''      
           &bands
 filband='Output_Si_band.dat',
/
 '''

# Eingabedatei wird erstellt   
nameIN_bands_mod = 'Input_Si_bands_mod.in'
file = open(nameIN_bands_mod,'w')
file.write(stringIN_bands_mod)
file.close()  

# Name der Ausgabedatei 
nameOUT_bands_mod = 'Ouput_Si_bands_mod.out'    

# Simulation wird durchgeführt
os.system('bands.x <'+nameIN_bands_mod+' >'+nameOUT_bands_mod)

#ermittlung der Fermienergie aus dem scf-Berechnung
with open(nameOUT_scf, 'r') as file:
        data = file.read()
        pos = data.find('highest occupied level (ev):     ') #Stelle finden an der die Fermienergie steht
    
    # Die Variable wird der Ausgabedatei entnommen
        fermienergy = floatFromString(data[pos+30:pos+50]) #herauslösen des Strings und konvertieren der Zahlen aus dem String in ein Float
        print(fermienergy)
       
stringIN_bands_plot = '''      
Output_Silizium_band.dat
-7 16
si.bands.xmgr
si.bands.ps
'''+str(fermienergy)+'''
1.0 '''+str(fermienergy)+'''
 '''

# Eingabedatei wird erstellt   
nameIN_bands_plot = 'Input_Si_bands_plot.in'
file = open(nameIN_bands_plot,'w')
file.write(stringIN_bands_plot)
file.close()  

# Name der Ausgabedatei 
nameOUT_bands_plot = 'Ouput_Si_bands_plot.out'    

# Simulation wird durchgeführt
os.system('plotband.x <'+nameIN_bands_plot+' >'+nameOUT_bands_plot)

fermiEnergy = 6.133 #eV

with open('Output_Silizium_band.dat', 'r') as f:
    d = f.readlines()
    i = 0
    data = np.array([])
    for line in d:
        if (i != 0):  # erste Zeile wird ignoriert            
            k = line.rstrip().split()
            if (i%2 == 1): # Zeile für Position im k-Raum
                arrayline = [] # neue Datenzeile wird angelegt
                for string in k:
                    arrayline.append(float(string))
            if (i%2 == 0): # Zeile für Energiebänder
                for string in k:
                    arrayline.append(float(string))
                arrayline = [np.array(arrayline)-fermiEnergy] # Verschiebung um Fermi-Energie
                if (len(data) == 0): # Erste Zeile im Datenarray wird angelegt
                    data = arrayline
                else: # falls nicht erstel Zeile werden weitere Zeilen angehängt
                    data = np.append(data,arrayline,axis =0)
        i+=1

f, ax = plt.subplots()
for i in range(3,11):
    plt.plot(data[:,i],'.', color='black',markersize=3)
    plt.plot(data[:,i], color='black')
plt.axhline(0,linestyle='--',color='green')
plt.xticks([])

#L 20
#gG 30
#X 10
#U 30
#gG 20

points = ['$L$','$\Gamma$','$X$','$U$','$\Gamma$']
width = [0,20,30,10,30]
ax.set_xlim([0,np.sum(width)])

ticks = []
pos = 0
for i in range(len(points)):
    pos = pos+width[i]
    plt.axvline(pos, color='black',linewidth=1.0)
    ticks.append(pos)
ax.set_xticks(ticks)
ax.set_xticklabels(points)
    
plt.ylabel('Energie (eV)')
plt.xlabel('Wellenvektor $k$')
plt.tick_params(axis='both', direction='in')

plt.savefig('Bandstruktur_Si.pdf')