#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Fortgeschrittenen Praktikum TU Berlin
# P60 Ramanspektroskopie an Halbleitern
# WS19/20 Gaebel & Krueger
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
import pickle
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
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
        ! Berechnungsart: SCF
     /
     &system
        ibrav = 4,
        ! Gitter-Struktur: hexagonal
        celldm(1) = 4.649,
        ! Gitterkonstante a in Bohr wird fetgelegt
        celldm(3)= 40,
        ! Der Abstand der Graphenlagen wir so groß gewählt, das faktisch nur eine Lage batrachtet wird
        nat = 2,
        ! Anzahl an Atomen pro Elementarzelle
        ntyp= 1,
        ! Anzahl an Atomsorten
        ecutwfc = 40,
        ! Energie, bis zu welcher die Funktionen entwickelt werden in Ry (CutOff)
        ecutrho  = 320,
	    ! maximale kinetische Energie der Elektronen (hier 8 x ecutwfc)
	nbnd = 8,
	! Anzahl der berechneten Bänder
     /
     &electrons
        mixing_beta = 0.7
        ! Gewichtung zwischen aufeinanderfolgenden Dichteverteilungen
     /
    ATOMIC_SPECIES
     C 12.011 C.pbe-rrkjus.UPF
     ! Element, Masse in a.m.u., Pseudopotential
    ATOMIC_POSITIONS crystal  
        C 0.0000000 0.0000000 0.000000  
        C 0.3333333 0.6666666 0.000000  
    K_POINTS crystal_b
     4
    	gG 40
        K 20
        M 30
        gG 10
    '''

# Eingabedatei wird erstellt   
nameIN_scf = 'Input_Graphen_scf.in'
file = open(nameIN_scf,'w')
file.write(stringIN_scf)
file.close()  

# Name der Ausgabedatei 
nameOUT_scf = 'Ouput_Graphen_scf.out'    

# Simulation wird durchgeführt
os.system('pw.x <'+nameIN_scf+' >'+nameOUT_scf)

stringIN_bands = '''      
 &control
        calculation = 'bands',
        ! Berechnungsart: SCF
     /
     &system
        ibrav = 4,
        ! Gitter-Struktur: hexagonal
        celldm(1) = 4.648726266579395,
        celldm(3)= 40,
        ! Gitterkonstante a in Bohr wird fetgelegt
        ! Der Abstand der Graphenlagen wir so groß gewählt, das faktisch nur eine Lage batrachtet wird
        nat = 2,
        ! Anzahl an Atomen pro Elementarzelle
        ntyp= 1,
        ! Anzahl an Atomsorten
        ecutwfc = 40,
        ! Energie, bis zu welcher die Funktionen entwickelt werden in Ry (CutOff)
        ecutrho  = 320,
	! maximale kinetische Energie der Elektronen (hier 8 x ecutwfc)
	nbnd = 8,
	! Anzahl der berechneten Bänder
     /
     &electrons
        mixing_beta = 0.7
        ! Gewichtung zwischen aufeinanderfolgenden Dichteverteilungen
     /
    ATOMIC_SPECIES
     C 12.011 C.pbe-rrkjus.UPF
     ! Element, Masse in a.m.u., Pseudopotential
    ATOMIC_POSITIONS crystal  
        C 0.0000000 0.0000000 0.000000  
        C 0.3333333 0.6666666 0.000000  
    K_POINTS crystal_b
     4
    	gG 40
        K 20
        M 30
        gG 10
 '''

# Eingabedatei wird erstellt   
nameIN_bands = 'Input_Graphen_bands.in'
file = open(nameIN_bands,'w')
file.write(stringIN_bands)
file.close()  

# Name der Ausgabedatei 
nameOUT_bands = 'Ouput_Graphen_bands.out'    

# Simulation wird durchgeführt
os.system('pw.x <'+nameIN_bands+' >'+nameOUT_bands)


stringIN_bands_mod = '''      
           &bands
 filband='Output_Graphen_band.dat',
/
 '''

# Eingabedatei wird erstellt   
nameIN_bands_mod = 'Input_Graphen_bands_mod.in'
file = open(nameIN_bands_mod,'w')
file.write(stringIN_bands_mod)
file.close()  

# Name der Ausgabedatei 
nameOUT_bands_mod = 'Ouput_Graphen_bands_mod.out'    

# Simulation wird durchgeführt
os.system('bands.x <'+nameIN_bands_mod+' >'+nameOUT_bands_mod)

#ermittlung der Fermienergie aus dem scf-Berechnung
with open('Ouput_Graphen_scf.out', 'r') as file:
        data = file.read()
        pos = data.find('highest occupied, lowest unoccupied level (ev):    ') #Stelle finden an der die Fermienergie steht
    
    # Die Variable wird der Ausgabedatei entnommen
        fermienergy = floatFromString(data[pos+50:pos+60]) #herauslösen des Strings und konvertieren der Zahlen aus dem String in ein Float
        print(fermienergy)
       
stringIN_bands_plot = '''      
Output_Graphen_band.dat
-7 16
Graphen.bands.xmgr
Graphen.bands.ps
'''+str(fermienergy)+'''
1.0 '''+str(fermienergy)+'''
 '''

# Eingabedatei wird erstellt   
nameIN_bands_plot = 'Input_Graphen_bands_plot.in'
file = open(nameIN_bands_plot,'w')
file.write(stringIN_bands_plot)
file.close()  

# Name der Ausgabedatei 
nameOUT_bands_plot = 'Ouput_Graphen_bands_plot.out'    

# Simulation wird durchgeführt
os.system('plotband.x <'+nameIN_bands_plot+' >'+nameOUT_bands_plot)


fermiEnergy =  -3.7939 #eV

with open('Output_Graphen_band.dat', 'r') as f:
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

points = ['$\Gamma$','$K$','$M$','$\Gamma$']
width = [0,40,20,30]
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

plt.savefig('Bandstruktur_Graphen.pdf')
pickle.dump(data,open('Data_Graphen_Band.p','wb'))
#os.system('rm *.in *.out *.wfc1 Graphen.bands.* *.dat *.gnu *.rap')