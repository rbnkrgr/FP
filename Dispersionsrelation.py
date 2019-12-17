#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Fortgeschrittenen Praktikum TU Berlin
# P60 Ramanspektroskopie an Halbleitern
# WS19/20 Gaebel & Krueger
"""
Berechnung der Dispersionsrelation
"""
import os
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

# Anzahl der verwendeten CPUs festlegen
CPU = 5

# Probe
probe = 'C_'

#---- Schritt 1  ----#
"""
     &control
        calculation = 'scf',
        ! Berechnungsart: SCF
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
        ecutwfc = 50,
        ! Energie, bis zu welcher die Funktionen entwickelt werden in Ry (CutOff)
     /
     &electrons
        ! mixing_beta = 0.7
        ! Gewichtung zwischen aufeinanderfolgenden Dichteverteilungen
	    conv_thr = 1e-10
	! Konvergenz-Schwellenwert für Energie 
     /
    ATOMIC_SPECIES
     Si  28.086 Si.rel-pbe-rrkj.UPF
     ! Element, Masse in a.m.u., Pseudopotential
    ATOMIC_POSITIONS alat
     Si 0.00 0.00 0.00
     Si 0.25 0.25 0.25   
     ! alat: Positionen werden in relativen Einheiten der Gitterkonstante angegeben
    K_POINTS automatic
         6  6  6  1  1  1
     ! automatic: Automatische Generierung der k-Punkte nach dem Monkhorst-System"""
     # Input Silizium     
     
input1 = """
       &control
        calculation = 'scf',
        ! Berechnungsart: SCF
     /
     &system
        ibrav = 4,
        ! Gitter-Struktur: hexagonal
        celldm(1) = 4.65,
        celldm(3)= 40,
        ! Gitterkonstante a in Bohr wird fetgelegt
        ! Der Abstand der Graphenlagen wir so groß gewählt, das faktisch nur eine Lage batrachtet wird
        nat = 2,
        ! Anzahl an Atomen pro Elementarzelle
        ntyp= 1,
        ! Anzahl an Atomsorten
        ecutwfc = 50,
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
        6  6  1  1  1  1"""
        # Input Graphene
     
# Eingabedatei wird erstellt   
nameIN = probe+'Input_Step1.in'
file = open(nameIN,'w')
file.write(input1)
file.close()  

# Name der Ausgabedatei 
nameOUT = probe+'Output_Step1.out'    

# Simulation wird durchgeführt
os.system('mpirun -np '+str(CPU)+' pw.x <'+nameIN+' | tee '+nameOUT)

#---- Schritt 2 ----#
"""
Phononen-Berechnung
&inputph
	tr2_ph = 1.0d-14,
	! Schwellenwert für Selbstkonsistenz  
	ldisp =.true.,
	! Phononendispersion wird bestimmt nach folgedendem Gitter
	nq1 = 4,
	nq2 = 4,
	nq3 = 4,
    ! Reziproke Gitterpunkte werden nach Monkhorst festgelegt
    ! da nur in 2D-Richtung untersucht, wird letzter Punkt auf 1 gesetzt
    fildyn = 'probe.dyn'
	! Outputfile: Datei mit der Dynamischen Matrix wird erstellt
/
0  0  0"""
# Input Silizium (Probe wird bei Bedarf ausgewählt)

input2 = """
Phononen-Berechnung
&inputph
	tr2_ph = 1.0d-14,
	! Schwellenwert für Selbstkonsistenz  
	ldisp =.true.,
	! Phononendispersion wird bestimmt nach folgedendem Gitter
	nq1 = 4,
	nq2 = 4,
	nq3 = 1,
    ! Reziproke Gitterpunkte werden nach Monkhorst festgelegt
    ! da nur in 2D-Richtung untersucht, wird letzter Punkt auf 1 gesetzt
    fildyn = 'probe.dyn'
	! Outputfile: Datei mit der Dynamischen Matrix wird erstellt
/
0  0  0"""
# Input Graphene

# Eingabedatei wird erstellt   
nameIN = probe+'Input_Step2.in'
file = open(nameIN,'w')
file.write(input2)
file.close()  

# Name der Ausgabedatei 
nameOUT = probe+'Output_Step2.out'    

# Simulation wird durchgeführt
os.system('mpirun -np '+str(CPU)+' ph.x <'+nameIN+' | tee '+nameOUT)

#---- Schritt 3 ----#
input3 = """
&input
	fildyn = 'probe.dyn',
	! Inputdatei
	zasr = 'simple',
	flfrc = 'ftout.fc',
	! Outputdatei
/
"""

# Eingabedatei wird erstellt   
nameIN = probe+'Input_Step3.in'
file = open(nameIN,'w')
file.write(input3)
file.close()  

# Name der Ausgabedatei 
nameOUT = probe+'Output_Step3.out'    

# Simulation wird durchgeführt
os.system('mpirun -np '+str(CPU)+' q2r.x <'+nameIN+' | tee '+nameOUT)

#---- Schritt 4 ----#
"""
&input
	asr = 'simple',
	flfrc = 'ftout.fc',
	! Input-Datei
	flfrq = 'Si.freq',
	! Ausgabe-Datei
	q_in_band_form = .true.,
	! nur Start- und End-Punkt im k-Raum müssen angegeben werden
	q_in_cryst_coord = .true.,
/
5
 gG 30 
 X  10
 K  20
 gG 30
 L  1
 """
# Input Silizium

input4 = """ 
&input
	asr = 'simple',
	flfrc = 'ftout.fc',
	! Input-Datei
	flfrq = 'C.freq',
	! Ausgabe-Datei
	q_in_band_form = .true.,
	! nur Start- und End-Punkt im k-Raum müssen angegeben werden
	q_in_cryst_coord = .true.,
/
4
 gG 30 
 M  20
 K  40
 gG 30
 """
# Input Graphene
 
# Eingabedatei wird erstellt   
nameIN = probe+'Input_Step4.in'
file = open(nameIN,'w')
file.write(input4)
file.close()  

# Name der Ausgabedatei 
nameOUT = probe+'Output_Step4.out'    

# Simulation wird durchgeführt
os.system('mpirun -np '+str(CPU)+' matdyn.x <'+nameIN+' | tee '+nameOUT)

#---- Schritt 5 ----#
"""
Si.freq
0 550
probe.phonon.bands.xmgr
probe.phonon.bands.ps
0
50 0
"""
# Input Silizium

input5 = """
C.freq
0 1700
probe.phonon.bands.xmgr
probe.phonon.bands.ps
0
50 0
"""
# Input Graphen

# Eingabedatei wird erstellt   
nameIN = probe+'Input_Step5.in'
file = open(nameIN,'w')
file.write(input5)
file.close()  

# Name der Ausgabedatei 
nameOUT = probe+'Output_Step5.out'    

# Simulation wird durchgeführt
os.system('plotband.x <'+nameIN+' >'+nameOUT)

#---- Schritt 6 ----#
# Plotten der Dispersionsrelation mit selbst erstelltem Programm

#-- Silizium --#
with open('Si.freq', 'r') as f:
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
                arrayline = [np.array(arrayline)] # Verschiebung um Fermi-Energie
                if (len(data) == 0): # Erste Zeile im Datenarray wird angelegt
                    data = arrayline
                else: # falls nicht erstel Zeile werden weitere Zeilen angehängt
                    data = np.append(data,arrayline,axis =0)
        i+=1

f, ax = plt.subplots()
for i in range(3,9):
    plt.plot(data[:,i],'.', color='black',markersize=3)
    plt.plot(data[:,i], color='black')
plt.xticks([])

#gG 30 
#X  10
#K  20
#gG 30
#L  1

points = ['$\Gamma$','$X$','$K$','$\Gamma$','$L$']
width = [0,30,10,20,30]
ax.set_xlim([0,np.sum(width)])

ticks = []
pos = 0
for i in range(len(points)):
    pos = pos+width[i]
    plt.axvline(pos, color='black',linewidth=1.0)
    ticks.append(pos)
ax.set_xticks(ticks)
ax.set_xticklabels(points)
    
plt.ylabel('Phononenfrequenz (1/cm)')
plt.xlabel('Wellenvektor $k$')
plt.tick_params(axis='both', direction='in')

#-- Graphen --#
with open('C.freq', 'r') as f:
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
                arrayline = [np.array(arrayline)] # Verschiebung um Fermi-Energie
                if (len(data) == 0): # Erste Zeile im Datenarray wird angelegt
                    data = arrayline
                else: # falls nicht erstel Zeile werden weitere Zeilen angehängt
                    data = np.append(data,arrayline,axis =0)
        i+=1

f, ax = plt.subplots()
for i in range(3,9):
    plt.plot(data[:,i],'.', color='black',markersize=3)
    plt.plot(data[:,i], color='black')
plt.xticks([])

#gG 30 
#M  20
#K  40
#gG 30

points = ['$\Gamma$','$M$','$K$','$\Gamma$']
width = [0,30,20,40]
ax.set_xlim([0,np.sum(width)])

ticks = []
pos = 0
for i in range(len(points)):
    pos = pos+width[i]
    plt.axvline(pos, color='black',linewidth=1.0)
    ticks.append(pos)
ax.set_xticks(ticks)
ax.set_xticklabels(points)
    
plt.ylabel('Phononenfrequenz (1/cm)')
plt.xlabel('Wellenvektor $k$')
plt.tick_params(axis='both', direction='in')