#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Fortgeschrittenen Praktikum TU Berlin
# P60 Ramanspektroskopie an Halbleitern
# WS19/20 Gaebel & Krueger
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from matplotlib import rc

rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

a0 = 0.529177210903 # Borscher Radius in Angström

CPU = 2 # Anzahl der verwendeten CPUs

def energyFit(a,e0,b,a0):
    return e0+b*(a-a0)**2

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
pathIN = ''
pathOUT = ''

#-- Graphen--#
energyList = []
# verschiedene Simulationen in Abhängigkeit des Parameters 'variable' werden durchgeführt
values= np.arange(2.36,2.56,0.02) # Hier CutOff
for variable in values:
    variableBohr = variable/a0# Umrechnung Angström zu Bohr
    print('Berechnung für a = '+str(variableBohr)+' Bohr')
    # Inputdatei wird generiert
    stringIN = '''  &control
        calculation = 'scf',
        ! Berechnungsart: SCF
     /
     &system
        ibrav = 4,
        ! Gitter-Struktur: hexagonal
        celldm(1) = '''+str(variableBohr)+''',
        celldm(3)= 40,
        ! Gitterkonstante a in Bohr wird fetgelegt
        ! Der Abstand der Graphenlagen wir so groß gewählt, das faktisch nur eine Lage batrachtet wird
        nat = 2,
        ! Anzahl an Atomen pro Elementarzelle
        ntyp= 1,
        ! Anzahl an Atomsorten
        ecutwfc = 40,
        ! Energie, bis zu welcher die Funktionen entwickelt werden in Ry (CutOff)
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
    K_POINTS automatic
        4  4  4  0  0  0
    
 '''

    # Eingabedatei wird erstellt   
    nameIN = pathIN+'Input_a'+str(variable)+'.in'
    file = open(nameIN,'w')
    file.write(stringIN)
    file.close()  
    
    # Name der Ausgabedatei 
    nameOUT = 'Output_a'+str(variable)+'.out'    

    # Simulation wird durchgeführt
    os.system('mpirun -np '+str(CPU)+' pw.x <'+nameIN+' | tee '+nameOUT)
    
    # Auswertung der Ausgabedatei
    fileOUT = open(pathOUT+nameOUT)
    
    with open(pathOUT+nameOUT, 'r') as file:
        data = file.read()
    pos = data.find('!    total energy              =     ')
    
    # Die Variable wird der Ausgabedatei entnommen
    if pos>0 :
        energy = floatFromString(data[pos+30:pos+50])
        print('Ergebnis: totale Energie = ' +str(energy)+' eV')
        energyList.append(energy)
    else:
        energyList.append(0)
    
    # Löschen der produizierten Dateien
    os.system('rm *.in *.out *.wfc1')   
    os.system('rm -r pwscf.save')
        
graphenList = energyList

#-- Diamant --#
energyList = []
# verschiedene Simulationen in Abhängigkeit des Parameters 'variable' werden durchgeführt
values= np.arange(3.46,3.66,0.02) # Hier CutOff
for variable in values:
    variableBohr = variable/a0# Umrechnung Angström zu Bohr
    print('Berechnung für a = '+str(variableBohr)+' Bohr')
    # Inputdatei wird generiert
    stringIN = '''  &control
        calculation = 'scf',
        ! Berechnungsart: SCF
     /
     &system
        ibrav = 2,
        ! Gitter-Struktur: fcc (diamant ist fcc mit 2 Basisatomen)
        celldm(1) = '''+str(variableBohr)+''',
        !celldm(2) = '''+str(variableBohr)+''',
        !celldm(3) = '''+str(variableBohr)+''',
        ! Gitterkonstante a in Bohr wird fetgelegt
        nat = 2,
        ! Anzahl an Atomen pro Elementarzelle
        ntyp= 1,
        ! Anzahl an Atomsorten
        ecutwfc = 180,
        ! Energie, bis zu welcher die Funktionen entwickelt werden in Ry (CutOff)
     /
     &electrons
        mixing_beta = 0.7
        ! Gewichtung zwischen aufeinanderfolgenden Dichteverteilungen
     /
    ATOMIC_SPECIES
     C 12.011 C.pbe-rrkjus.UPF
     ! Element, Masse in a.m.u., Pseudopotential
      ATOMIC_POSITIONS alat
     C 0.00 0.00 0.00
     C 0.25 0.25 0.25
     ! alat: Positionen werden in relativen Einheiten der Gitterkonstante angegeben
    K_POINTS automatic
        4  4  4  0  0  0
    
 '''

    # Eingabedatei wird erstellt   
    nameIN = pathIN+'Input_a'+str(variable)+'.in'
    file = open(nameIN,'w')
    file.write(stringIN)
    file.close()  
    
    # Name der Ausgabedatei 
    nameOUT = 'Output_a'+str(variable)+'.out'    

    # Simulation wird durchgeführt
    os.system('mpirun -np '+str(CPU)+' pw.x <'+nameIN+' | tee '+nameOUT)
    
    # Auswertung der Ausgabedatei
    fileOUT = open(pathOUT+nameOUT)
    
    with open(pathOUT+nameOUT, 'r') as file:
        data = file.read()
    pos = data.find('!    total energy              =     ')
    
    # Die Variable wird der Ausgabedatei entnommen
    if pos>0 :
        energy = floatFromString(data[pos+30:pos+50])
        print('Ergebnis: totale Energie = ' +str(energy)+' eV')
        energyList.append(energy)
    else:
        energyList.append(0)
   
    # Löschen der produizierten Dateien
    os.system('rm *.in *.out *.wfc1')   
    os.system('rm -r pwscf.save')
        
diamantList = energyList

#pickle.dump((values,energyList), open('GitterC.p','wb'))
#values,energyList = pickle.load( open('GitterC.p','rb'))

fig1, ax = plt.subplots(1, 2, sharey ='row')
plt.subplots_adjust(wspace=0,hspace=0) # kein Platz zwischen einzelnen Subplots

ax[0].plot(np.arange(2.36,2.56,0.02),graphenList,'.',color='black')
params, params_covariance = optimize.curve_fit(energyFit,np.arange(2.36,2.56,0.02),graphenList)
ax[0].plot(np.arange(2.345,2.57,0.001),energyFit(np.arange(2.345,2.57,0.001),params[0],params[1],params[2]),color='green',label='Graphen')
ax[0].set_ylabel('totale Energie (eV)')
ax[0].tick_params(axis='both', direction='in')
ax[0].set_title('Graphen')
ax[0].set_xlabel('Gitterkonstante (Å)')

ax[1].plot(values,energyList,'.',color='black')
params, params_covariance = optimize.curve_fit(energyFit,values,energyList)
ax[1].plot(np.arange(3.45,3.67,0.001),energyFit(np.arange(3.45,3.67,0.001),params[0],params[1],params[2]),color='blue',label='Diamant')
ax[1].tick_params(axis='both', direction='in')
ax[1].set_title('Diamant')
ax[1].set_xlabel('Gitterkonstante (Å)')
plt.savefig('CGitterKonst.pdf')

#-- Silizium --#
energyList = []
# verschiedene Simulationen in Abhängigkeit des Parameters 'variable' werden durchgeführt
values= np.arange(5.35,5.55,0.02) # Hier CutOff
for variable in values:
    variableBohr = variable/a0# Umrechnung Angström zu Bohr
    print(variableBohr)
    # Inputdatei wird generiert
    stringIN = '''  &control
        calculation = 'scf',
        ! Berechnungsart: SCF
     /
     &system
        ibrav = 2,
        ! Gitter-Struktur: fcc (diamant ist fcc mit 2 Basisatomen)
        celldm(1) = '''+str(variableBohr)+''',
        !celldm(2) = '''+str(variableBohr)+''',
        !celldm(3) = '''+str(variableBohr)+''',
        ! Gitterkonstante a in Bohr wird fetgelegt
        nat = 2,
        ! Anzahl an Atomen pro Elementarzelle
        ntyp= 1,
        ! Anzahl an Atomsorten
        ecutwfc = 180,
        ! Energie, bis zu welcher die Funktionen entwickelt werden in Ry (CutOff)
     /
     &electrons
        mixing_beta = 0.7
        ! Gewichtung zwischen aufeinanderfolgenden Dichteverteilungen
     /
    ATOMIC_SPECIES
     Si  28.086 Si.pbe-rrkj.UPF
     ! Element, Masse in a.m.u., Pseudopotential
      ATOMIC_POSITIONS alat
     Si 0.00 0.00 0.00
     Si 0.25 0.25 0.25
     ! alat: Positionen werden in relativen Einheiten der Gitterkonstante angegeben
     !   K_POINTS 
     !1
     ! 0.0  0.0  0.0   1.0
      K_POINTS automatic
         8  8  8  0  0  0
 '''

    # Eingabedatei wird erstellt   
    nameIN = pathIN+'Input_a'+str(variable)+'.in'
    file = open(nameIN,'w')
    file.write(stringIN)
    file.close()  
    
    # Name der Ausgabedatei 
    nameOUT = 'Output_a'+str(variable)+'.out'    

    # Simulation wird durchgeführt
    os.system('mpirun -np '+str(CPU)+' pw.x <'+nameIN+' | tee '+nameOUT)
    
    # Auswertung der Ausgabedatei
    fileOUT = open(pathOUT+nameOUT)
    
    with open(pathOUT+nameOUT, 'r') as file:
        data = file.read()
    pos = data.find('!    total energy              =     ')
    
    # Die Variable wird der Ausgabedatei entnommen
    if pos>0 :
        energy = floatFromString(data[pos+30:pos+50])
        energyList.append(energy)
    else:
        energyList.append(0)

    # Löschen der produizerten Dateien
    #os.system('rm *.in *.out')   
    os.system('rm -r pwscf.save')
        
plt.plot(values,energyList,'.',color='black')
params, params_covariance = optimize.curve_fit(energyFit,values,energyList)
plt.plot(np.arange(5.34,5.56,0.001),energyFit(np.arange(5.34,5.56,0.001),params[0],params[1],params[2]),color='green')
plt.tick_params(axis='both', direction='in')
plt.xlabel('Gitterkonstante (Å)')
plt.ylabel('totale Energie (eV)')

plt.savefig('SiGitterKonst.pdf')