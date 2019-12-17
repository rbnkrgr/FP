#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Fortgeschrittenen Praktikum TU Berlin
# P60 Ramanspektroskopie an Halbleitern
# WS19/20 Gaebel & Krueger
import os
import time
import matplotlib.pyplot as plt
import numpy as np
import pickle

# Borscher radius in Angström
a0 = 0.529177210903 

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

# Einführung der Ausgabelisten
energyList = []
timeList   = []
timeList1  = []
RamList    = []
forceList  = []

# verschiedene Simulationen in Abhängigkeit des Parameters 'variable' werden durchgeführt 
values1 = np.arange(10,300,1) # Hier CutOff-Liste erstellen
for variable in values1:
    
    # Umrechnung Angström zu Bohr
    variableBohr = 5.43/a0
    
    print(variable)
    
    # Inputdatei wird generiert
    stringIN = '''  &control
        calculation = 'scf',
        ! Berechnungsart: SCF
        tprnfor = .TRUE.
        ! Die Kraftkomponenten werden berechnet
     /
     &system
        ibrav = 2,
        ! Gitter-Struktur: fcc (diamant ist fcc mit 2 Basisatomen)
        celldm(1) = '''+str(variableBohr)+''',
        ! Gitterkonstante a in Bohr wird fetgelegt
        nat = 2,
        ! Anzahl an Atomen pro Elementarzelle
        ntyp= 1,
        ! Anzahl an Atomsorten
        ecutwfc = '''+str(variable)+''',
        ! Energie, bis zu welcher die Funktionen entwickelt werden in Ry (CutOff)
     /
     &electrons
        mixing_beta = 0.7
        ! Gewichtung zwischen aufeinanderfolgenden Dichteverteilungen
      /
     ATOMIC_SPECIES
     Si  28.086 Si.pz-vbc.UPF
     ! Element, Masse in a.m.u., Pseudopotential
     ATOMIC_POSITIONS alat
     Si 0.00 0.00 0.00
     Si 0.25 0.25 0.25   
     ! alat: Positionen werden in relativen Einheiten der Gitterkonstante angegeben
    ! K_POINTS automatic
     !    4  4  4  0  0  0
      K_POINTS
       1
        0.0  0.0  0.0   1.0
 '''

    # Eingabedatei wird erstellt   
    nameIN = pathIN+'Input_CO'+str(variable)+'.in'
    file = open(nameIN,'w')
    file.write(stringIN)
    file.close()  
    
    # Name der Ausgabedatei 
    nameOUT = 'Output_CO'+str(variable)+'.out'    

    #zeit start
    t1 = time.time()
    # Simulation wird durchgeführt
    os.system('mpirun -np 4 '+pathQE+' <'+pathIN+nameIN+' >'+pathOUT+nameOUT)
    #zeit ende
    t2 = time.time()
    
    # Bestimmung der Ausführungsdauer des Programms
    dt = t2-t1
    print(dt)
    timeList.append(dt)
    
    # Auswertung der Ausgabedatei
    fileOUT = open(pathOUT+nameOUT)
    
    with open(pathOUT+nameOUT, 'r') as file:
        data = file.read()
    pos = data.find('!    total energy              =     ')
    timepos = data.find('     PWSCF        :')
    Rampos = data.find('RAM per process > ')
    forcepos = data.find('Total force =')
    
    # Die Variable wird der Ausgabedatei entnommen
    time1 = floatFromString(data[timepos+17:timepos+40])
    timeList1.append(time1)
    Ram = floatFromString(data[Rampos+16:Rampos+29])
    RamList.append(Ram)
    force = floatFromString(data[forcepos+16:forcepos+29])
    forceList.append(force)
    if pos>0 :
        energy = floatFromString(data[pos+30:pos+50])
        energyList.append(energy)
    else:
        energyList.append(0)    

    pickle.dump((energyList,values1,timeList1,timeList,RamList,forceList),open('Save-CutOff_1_800.p','wb'))

    #erneutes erstellen des Inputs mit größerem Abstand der Cutoff-Werte
values2 = np.arange(300,800,5) # Hier CutOff-Liste erstellen
for variable in values2:
    
    # Umrechnung Angström zu Bohr
    variableBohr = 5.43/a0
    
    print(variable)
    
    # Inputdatei wird generiert
    stringIN = '''  &control
        calculation = 'scf',
        ! Berechnungsart: SCF
        tprnfor = .TRUE.
        ! Die Kraftkomponenten werden berechnet
     /
     &system
        ibrav = 2,
        ! Gitter-Struktur: fcc (diamant ist fcc mit 2 Basisatomen)
        celldm(1) = '''+str(variableBohr)+''',
        ! Gitterkonstante a in Bohr wird fetgelegt
        nat = 2,
        ! Anzahl an Atomen pro Elementarzelle
        ntyp= 1,
        ! Anzahl an Atomsorten
        ecutwfc = '''+str(variable)+''',
        ! Energie, bis zu welcher die Funktionen entwickelt werden in Ry (CutOff)
     /
     &electrons
        mixing_beta = 0.7
        ! Gewichtung zwischen aufeinanderfolgenden Dichteverteilungen
      /
     ATOMIC_SPECIES
     Si  28.086 Si.pz-vbc.UPF
     ! Element, Masse in a.m.u., Pseudopotential
     ATOMIC_POSITIONS alat
     Si 0.00 0.00 0.00
     Si 0.25 0.25 0.25   
     ! alat: Positionen werden in relativen Einheiten der Gitterkonstante angegeben
    ! K_POINTS automatic
     !    4  4  4  0  0  0
      K_POINTS
       1
        0.0  0.0  0.0   1.0
 '''

    # Eingabedatei wird erstellt   
    nameIN = pathIN+'Input_CO'+str(variable)+'.in'
    file = open(nameIN,'w')
    file.write(stringIN)
    file.close()  
    
    # Name der Ausgabedatei 
    nameOUT = 'Output_CO'+str(variable)+'.out'    
    
    t1 = time.time()
    # Simulation wird durchgeführt
    os.system('mpirun -np 4 '+pathQE+' <'+pathIN+nameIN+' >'+pathOUT+nameOUT)
    t2 = time.time()
    
    # Bestimmung der Ausführungsdauer des Programms
    dt = t2-t1
    print(dt)
    timeList.append(dt)
    
    # Auswertung der Ausgabedatei
    fileOUT = open(pathOUT+nameOUT)
    
    with open(pathOUT+nameOUT, 'r') as file:
        data = file.read()
    pos = data.find('!    total energy              =     ')
    timepos = data.find('     PWSCF        :')
    Rampos = data.find('RAM per process > ')
    forcepos = data.find('Total force =')
    
    # Die Variable wird der Ausgabedatei entnommen
    time1 = floatFromString(data[timepos+17:timepos+40])
    timeList1.append(time1)
    Ram = floatFromString(data[Rampos+16:Rampos+29])
    RamList.append(Ram)
    force = floatFromString(data[forcepos+16:forcepos+29])
    forceList.append(force)
    if pos>0 :
        energy = floatFromString(data[pos+30:pos+50])
        energyList.append(energy)
    else:
        energyList.append(0)    
    values = np.append(values1,values2)

pickle.dump((energyList,values,timeList1,timeList,RamList,forceList),open('Save-CutOff_1_800.p','wb'))
#energyList,values,timeList1,timeList,RamList,forceList = pickle.load(open('Save-CutOff_300_800.p','rb'))

energyPerAtom = np.array(energyList)/2

fig, ax1 = plt.subplots()
ax1.plot(values,energyPerAtom,'.')
ax2 = ax1.twinx()
ableitung = np.gradient(energyList)
ax2.plot(values,10**3*ableitung)

# Plotten
fig1, ax1 = plt.subplots(2, sharey='row', sharex ='col') # Gemeinsame Achsen    
plt.subplots_adjust(wspace=0,hspace=0) # Kein Platz zwischen einzelnen Subplots

lns1 = ax1[0].plot(values,timeList1,'s',color='red',label='Zeit', markersize = 1) # Plotten der Messdaten
ax12 = ax1[0].twinx()
lns2 = ax12.plot(values,RamList,'v',color='green',label='RAM', markersize = 1) # Plotten der Messdaten
lns = lns1+lns2
labs = [l.get_label() for l in lns]
ax1[0].legend(lns, labs, loc=0)
lns3 = ax1[1].plot(values,energyList,'^',color='blue',label='Energie', markersize = 1) # Plotten der Messdaten  
ax22 = ax1[1].twinx()
lns4 = ax22.plot(values,forceList,'o',color='black',label='Kraft', markersize = 1)    
lns = lns3+lns4
labs = [l.get_label() for l in lns]
ax1[1].legend(lns, labs, loc=0)

ax1[0].tick_params(axis='both', direction='in')
ax1[1].tick_params(axis='both', direction='in')
ax12.tick_params(axis='both', direction='in')
ax22.tick_params(axis='both', direction='in')             

ax12.set_ylabel('Arbeitsspeicher (MB)')
ax1[0].set_ylabel('CPU-Zeit (s)')
ax1[1].set_ylabel('totale Energie (eV)')
ax1[1].set_xlabel('mesh cutoff Wert (Ry)')
ax22.set_ylabel('max. Kraftkomponente')

plt.savefig('Uebersicht.pdf', bbox_inches="tight")

# Löschen der erzeugten Dateien
os.system('rm -r *.in *.out *.wfc1 *.xml *.save *.mix1')