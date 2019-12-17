#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Fortgeschrittenen Praktikum TU Berlin
# P60 Ramanspektroskopie an Halbleitern
# WS19/20 Gaebel & Krueger
import os

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

#--- Silizium ---#

# Inputdatei wird generiert
stringIN_scf = '''      
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
     ! automatic: Automatische Generierung der k-Punkte nach dem Monkhorst-Packdfs
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

# Inputdatei wird generiert
stringIN_scf = '''      
 Phononen-Berechnung am Gammapunkt
&inputph
	tr2_ph = 1.0d-14,
	! Schwellenwert für Selbstkonsistenz  
	epsil =.true.,
	! Phononendispersion wird bestimmt nach folgedendem Gitter
	! Reziproke Gitterpunkte werden nach Monkhorst festgelegt
	fildyn = 'si.dyn'
	! Outputfile: Datei mit einer Dynamischen Matrix wird erstellt
/
0  0  0
 '''

# Eingabedatei wird erstellt   
nameIN_scf = 'Input_Si_scf.in'
file = open(nameIN_scf,'w')
file.write(stringIN_scf)
file.close()  

# Name der Ausgabedatei 
nameOUT_scf = 'Ouput_Si_scf.out'    

# Simulation wird durchgeführt
os.system('ph.x <'+nameIN_scf+' >'+nameOUT_scf) 

peaks_thz_Si = []
peaks_cm_Si = []
print('Peaks für Silizium:')
for i in range (1,7):
    with open(nameOUT_scf, 'r') as file:
        data = file.read()
        pos = data.find('     freq (    '+str(i)+') =') #Stelle finden an der die Fermienergie steht
    
    # Die Variable wird der Ausgabedatei entnommen
        freq_thz = floatFromString(data[pos+20:pos+35])
        freq_cm = floatFromString(data[pos+40:pos+58])#herauslösen des Strings und konvertieren der Zahlen aus dem String in ein Float
        print(freq_thz, freq_cm)
        peaks_thz_Si.append(freq_thz)
        peaks_cm_Si.append(freq_cm)

#--- Graphen ---#

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
     /
     &electrons
        ! mixing_beta = 0.7
        ! Gewichtung zwischen aufeinanderfolgenden Dichteverteilungen
	    conv_thr = 1e-10
	     ! Konvergenz-Schwellenwert für Energie 
     /
     ATOMIC_SPECIES
     C 12.011 C.pbe-rrkjus.UPF
     ! Element, Masse in a.m.u., Pseudopotential
    ATOMIC_POSITIONS crystal  
        C 0.0000000 0.0000000 0.000000  
        C 0.3333333 0.6666666 0.000000  
    K_POINTS automatic
         6  6  6  1  1  1
     ! automatic: Automatische Generierung der k-Punkte nach dem Monkhorst-Packdfs
 '''

# Eingabedatei wird erstellt   
nameIN_scf = 'Input_Graphen_scf.in'
file = open(nameIN_scf,'w')
file.write(stringIN_scf)
file.close()  

# Name der Ausgabedatei 
nameOUT_scf = 'Ouput_Graphen_scf.out'    

# Simulation wird durchgeführt
os.system('pw.x <'+nameIN_scf+' > '+nameOUT_scf)

# Inputdatei wird generiert
stringIN_scf = '''      
 Phononen-Berechnung am Gammapunkt
&inputph
	tr2_ph = 1.0d-14,
	! Schwellenwert für Selbstkonsistenz  
	epsil =.true.,
	! Phononendispersion wird bestimmt nach folgedendem Gitter
	! Reszproke Gitterpunkte werden nach Monkhorst festgelegt
	fildyn = 'si.dyn'
	! Outputfile: Datei mit einer Dynamischen Matrix wird erstellt
/
0  0  0
 '''

# Eingabedatei wird erstellt   
nameIN_scf = 'Input_Graphen_scf.in'
file = open(nameIN_scf,'w')
file.write(stringIN_scf)
file.close()  

# Name der Ausgabedatei 
nameOUT_scf = 'Ouput_Graphen_scf.out'    

# Simulation wird durchgeführt
os.system('ph.x <'+nameIN_scf+' >'+nameOUT_scf) 

peaks_thz_Graphen = []
peaks_cm_Graphen = []
print('Peaks für Graphen:')
for i in range (1,7):
    with open(nameOUT_scf, 'r') as file:
        data = file.read()
        pos = data.find('     freq (    '+str(i)+') =') #Stelle finden an der die Fermienergie steht
    
    # Die Variable wird der Ausgabedatei entnommen
        freq_thz = floatFromString(data[pos+20:pos+35])
        freq_cm = floatFromString(data[pos+40:pos+58])#herauslösen des Strings und konvertieren der Zahlen aus dem String in ein Float
        print(freq_thz, freq_cm)
        peaks_thz_Graphen.append(freq_thz)
        peaks_cm_Graphen.append(freq_cm)
        
#pickle.dump((peaks_cm_Graphen, peaks_thz_Graphen, peaks_cm_Si, peaks_thz_Si),open('Peak_point.p','wb'))
#peaks_cm_Graphen, peaks_thz_Graphen, peaks_cm_Si, peaks_thz_Si = pickle.load(open('Peak_point.p','rb'))