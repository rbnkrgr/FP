# -*- coding: utf-8 -*-
# Fortgeschrittenen Praktikum TU Berlin
# P60 Ramanspektroskopie an Halbleitern
# WS19/20 Simon Gaebel & Robin Krüger

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import odr
import os
import pandas as pd

# Festlegen der Schriftarten für matplotlib
plt.rcParams['font.sans-serif'] = "Dubai"
plt.rcParams['font.family'] = "sans-serif"

def Lorentzian(x, amp, cen, wid):
    '''Definition der Lorentz-Funktion zum Fitten'''
    return (amp*wid**2/((x-cen)**2+wid**2))

def LorentzianOff(x, amp, cen, wid, off):
    '''Definition der Lorentz-Funktion zum Fitten mit Offset'''
    return (amp*wid**2/((x-cen)**2+wid**2)+off)

def sinSquare(beta,x):
    '''Definition der Funktion zum Fitten der
    Messungen mit paralleler Polarisation -> Darstellung für ODR'''
    amp,phi = beta
    return amp*(np.sin(2*(x-phi))**2)

def cosSquare(beta, x):
    '''Definition der Funktion zum Fitten der
    Messungen mit paralleler Polarisation -> Darstellung für ODR'''
    amp, phi = beta
    return amp*(np.cos(2*(x-phi))**2)

def ramanstyle():
    '''Festlegung von Formatierung und Beschriftung der Raman-Plots'''
    plt.tick_params(axis='both', direction='in')
    plt.xlabel('Ramanverschiebung (1/cm)')
    plt.ylabel('Intensität (arb. u.)')

def giveArea(data,min,max):
    '''Ermittlung des Bereiches eines Arrays
    -> Spezifikation über Intervall [min,max] in der ersten Spalte'''
    area = np.argwhere(data[:,0]>=min)
    min_index = int(area[0])
    area = np.argwhere(data[:,0]>=max)
    max_index = int(area[0])
    return data[min_index:max_index,:]

###---------- Kalibriermessung Si-111 ----------###
# Datenimport
daten = np.loadtxt('kalibrierung_allg.txt')

# Festlegen der bertachteten Grenzen
min = 400
max = 600
ausschnittDarst = giveArea(daten,200,1100)
ausschnittFit = giveArea(daten,min,max)

# Fitten des Peaks und Plotten
params, params_covariance = optimize.curve_fit(Lorentzian,ausschnittFit[:,0],ausschnittFit[:,1],p0=[1,517,1])
plt.plot(np.arange(min,max,0.2),Lorentzian(np.arange(min,max,0.2),params[0],params[1],params[2]),color='red',linewidth=0.9)
plt.plot(ausschnittDarst[:,0],ausschnittDarst[:,1],'.',markersize=2,color='black')
ramanstyle()
plt.savefig('Kalibrierung3.pdf')
plt.show()

# Ausgabe der Parameter
cen_err = np.sqrt(params_covariance[1,1])
print('Der Peak hat die Position: ('+str(params[1])+'+-'+str(cen_err)+') 1/cm')
verschiebung = 520 - params[1]
print('Somit werden die Folgenden Spektren verschoben um: '+str(verschiebung)+' 1/cm')

###---------- Orientierung der Probe Si-100 ----------###
# Pfad zu Messdaten und Ergebnisspeicherung
mPath = 'Silizium'
ePath = 'Ergebnis'

# Betrachteter Bereich der Verschiebung
min = 475
max = 560

# DataFrame-Objekt (pandas) zur Speicherung der Messwerte und Fit-Parameter 
d = {'polWinkel':[],'winkel':[],'datenAusschnitt':[],'fPamp':[],'fPamp_var':[],'fPcen':[],'fPcen_var':[],'fPwid':[],'fPwid_var':[]}
theObject = pd.DataFrame(data=d)

for filename in os.listdir(mPath):
    '''Einlesen der Messdateien und Auslesen der Paramter
    Dateinamenform: Si100_A_000_00.txt
    000 -> Probenwinkel, 00 -> eingestellter Winkel am Polarisator'''
    daten = np.loadtxt(mPath+'/'+filename)
    winkel = int(filename[8:11])
    polWinkel = int(filename[12:14])
    
    # Wahl des betrachteten Ausschnitts und Fit
    ausschnitt = giveArea(daten,475,560)
    params, params_covariance = optimize.curve_fit(Lorentzian,ausschnitt[:,0],ausschnitt[:,1],p0=[1,517,1])
    
    # Plotten der einzelnen Fits und abspeichern der Bilder
    title = 'Silizium-100, Winkel: '+str(winkel)+' Grad, Polarisation: '+str(polWinkel*2)+' Grad'
    plt.plot(np.arange(min,max,0.2)+verschiebung,Lorentzian(np.arange(min,max,0.2),params[0],params[1],params[2]),color='red',linewidth=0.9)
    plt.plot(ausschnitt[:,0]+verschiebung,ausschnitt[:,1],'.',markersize=2,color='black')
    ramanstyle()
    plt.savefig(ePath+'/'+filename[0:14]+'.pdf')
    plt.show()

    # Kombinierung der Mess- und Fit-Werte in einem DataFrame-Object (pandas)
    d = {'polWinkel':polWinkel,'winkel':winkel,'datenAusschnitt':[ausschnitt],'fPamp':params[0],'fPamp_var':params_covariance[0,0],'fPcen':params[1],'fPcen_var':params_covariance[1,1],'fPwid':params[2],'fPwid_var':params_covariance[2,2]}
    appendObject = pd.DataFrame(data=d)
    
    # Anängen der Daten an DataFrame-Objekt zur Datensammlung
    theObject = theObject.append(appendObject)
    
# Die gesammelten Daten werden neu indiziert
theObject = theObject.set_index(np.arange(theObject.shape[0]))

# Die Zeilen werden erst nach Polariastion und dann nach Winkel sortiert
sortiert = theObject.sort_values(by=['polWinkel','winkel'])

# Trennung der Daten für beide Polarisationsrichtungen
sortiertPol0  = sortiert.take(np.arange(19),axis=0)
sortiertPol45 = sortiert.take(np.arange(19,38),axis=0)

# Daten zum Fitten und Plotten auslesen
datenauswahl = sortiertPol0.take([1,3,4],axis=1)
datenarray0 = datenauswahl.to_numpy()
datenauswahl = sortiertPol45.take([1,3,4],axis=1)
datenarray45 = datenauswahl.to_numpy()

# Plotten
winkelRad = np.deg2rad(datenarray0[:,0])

# Fitten mit Orthogonal Distance Regression
# Daten aus paralleler Polarisation
Model0 = odr.Model(sinSquare)
Data0 = odr.RealData(x=winkelRad,y=datenarray0[:,1],sx=np.deg2rad(1),sy=np.sqrt(np.abs(datenarray0[:,2])))
Fit0 = odr.ODR(Data0,Model0,beta0=[1,0])
Out0 = Fit0.run() 
Out0.pprint()  # Ergebnisse des Fits in Out0.beta und Fehler in Out0.sd_beta

# Daten aus senkrechter Polarisation
Model45 = odr.Model(cosSquare)
Data45 = odr.RealData(x=winkelRad,y=datenarray45[:,1],sx=np.deg2rad(1),sy=np.sqrt(np.abs(datenarray45[:,2])))
Fit45 = odr.ODR(Data45,Model45,beta0=[1,0])
Out45 = Fit45.run() 
Out45.pprint()

# Plotten
fig1, ax = plt.subplots(2, 1, sharex ='col')
plt.subplots_adjust(wspace=0,hspace=0) # kein Platz zwischen einzelnen Subplots

# parallele Polarisation
ax[0].errorbar(winkelRad,datenarray0[:,1],xerr=np.deg2rad(1),yerr=np.sqrt(np.abs(datenarray0[:,2])),ls=' ',color='black')
ax[0].plot(np.arange(0,1.1*np.pi,0.05),sinSquare(Out0.beta,np.arange(0,1.1*np.pi,0.05)),color='g',linewidth=0.9)
ax[0].set_ylabel('parallele Polarisation')
ax[0].tick_params(axis='both', direction='in')
ax[0].set_yticks([]) # Unterdrücken der Ticks

# senkrechte Polarisation
ax[1].errorbar(winkelRad,datenarray45[:,1],xerr=np.deg2rad(1),yerr=np.sqrt(np.abs(datenarray45[:,2])),ls=' ',color='black')
ax[1].plot(np.arange(0,1.1*np.pi,0.05),cosSquare(Out45.beta,np.arange(0,1.1*np.pi,0.05)),color='g',linewidth=0.9)
ax[1].set_xlabel('Ausrichtungswinkel der Probe (rad)')
ax[1].set_ylabel('senkrechte Polarisation')
ax[1].tick_params(axis='both', direction='in')
ax[1].set_yticks([])
ax[1].text(-0.5,400,'Amplitude der Lorentzfunktion (arb. u.)',rotation=90)

plt.savefig(ePath+'/Orientierung.pdf')
plt.show()

# Berechnung des Vershiebungswinkels
phi1 = np.rad2deg(Out0.beta[1])
phi1_err = np.rad2deg(Out0.sd_beta[1])
phi2 = np.rad2deg(Out45.beta[1]%(np.pi))
phi2_err = np.rad2deg(Out45.sd_beta[1])
phiMittel = (phi1+phi2)/2
phiMittel_err = 0.5*np.sqrt(phi1_err**2+phi2_err**2)
print('Phi bei paralleler Polarisation:'+str(phi1)+'+-'+str(phi1_err))
print('Phi bei senkrechter Polarisation:'+str(phi2)+'+-'+str(phi2_err))
print('Phi im Mittel:'+str(phiMittel)+'+-'+str(phiMittel_err))

###---------- Auswertung Graphit ----------###
# Festlegen des Quellpfades
mPath = 'Graphit'

# DataFrame-Objekt (pandas) zur Speicherung der Messwerte -> Wird hier verwendet, damit alle Daten gesammlt in Python vorliegen
d = {'Name':[],'Daten':[],'K1':[],'K1cov':[],'K2':[],'K2cov':[],'K3':[],'K3cov':[]}
graphitObject = pd.DataFrame(data=d)

for filename in os.listdir(mPath):
    daten = np.loadtxt(mPath+'/'+filename)
    
    # Anängen der Daten an DataFrame-Objekt zur Datensammlung
    d = {'Name':filename[0:-4],'Daten':[daten]}
    appendObject = pd.DataFrame(data=d)
    graphitObject = graphitObject.append(appendObject)
    
# Die gesammelten Daten werden neu indiziert
graphitObject = graphitObject.set_index(np.arange(graphitObject.shape[0]))

# Betrachteter Bereich
min = 1100
max = 2900

# Fit-Bereiche für einzelne Peaks
minK1 = 1250
maxK1 = 1450
minK2 = 1400
maxK2 = 2250
minK3 = 2250
maxK3 = 2800

# Indizes der Proben festlegen P0 (0+5), P1 (2+7), P2 (3+8), P3 (4+9) (nach diesen sind sie in graphitObject aufgeführt)
Messungen = [[0,5],[2,7],[3,8],[4,9]]
fig1, ax1 = plt.subplots(4, 2, sharey='row', sharex ='col') # Gemeinsame Achsen
plt.subplots_adjust(wspace=0,hspace=0) # Kein Platz zwischen einzelnen Subplots

for Zeile in range(4):
    for Spalte in range(2):
        i = Messungen[Zeile][Spalte]   
        data_im = graphitObject.at[i,'Daten']
        datav = np.transpose([data_im[:,0]+verschiebung,data_im[:,1]]) # Verschieben der Spektren um Wert aus Kalibriermessung
        datav2 = giveArea(datav,min,max) # Betrachtung eines Ausschnittes der Daten
        ax1[Zeile,Spalte].plot(datav2[:,0],datav2[:,1],'.',markersize=2,color='black') # Plotten der Messdaten
       
        # Anlegen der einzelnen Fits
        # K1
        ausschnittK1 = giveArea(datav2,minK1,maxK1) # Auschnitt zum Fitten
        paramsK1, params_covarianceK1 = optimize.curve_fit(LorentzianOff,ausschnittK1[:,0],ausschnittK1[:,1],p0=[1,1330,1,4]) # Fitten
        ax1[Zeile,Spalte].plot(datav2[:,0],LorentzianOff(datav2[:,0],paramsK1[0],paramsK1[1],paramsK1[2],paramsK1[3]),color='g',linewidth=1) # Plotten
        graphitObject.at[i,'K1'] = paramsK1[0] # Amplitude und Varianz werden zur Berechnung der Depolarisationsrate gespeichert.
        graphitObject.at[i,'K1cov'] = params_covarianceK1[0,0]
        
        # K2
        ausschnittK2 = giveArea(datav2,minK2,maxK2)
        paramsK2, params_covarianceK2 = optimize.curve_fit(LorentzianOff,ausschnittK2[:,0],ausschnittK2[:,1],p0=[1,1575,1,4])
        ax1[Zeile,Spalte].plot(datav2[:,0],LorentzianOff(datav2[:,0],paramsK2[0],paramsK2[1],paramsK2[2],paramsK2[3]),color='g',linewidth=1)
        graphitObject.at[i,'K2'] = paramsK2[0]
        graphitObject.at[i,'K2cov'] = params_covarianceK2[0,0]
        
        # K3
        ausschnittK3 = giveArea(datav2,minK3,maxK3)
        paramsK3, params_covarianceK3 = optimize.curve_fit(LorentzianOff,ausschnittK3[:,0],ausschnittK3[:,1],p0=[1,2650,1,4])
        ax1[Zeile,Spalte].plot(datav2[:,0],LorentzianOff(datav2[:,0],paramsK3[0],paramsK3[1],paramsK3[2],paramsK3[3]),color='g',linewidth=1)
        graphitObject.at[i,'K3'] = paramsK3[0]
        graphitObject.at[i,'K3cov'] = params_covarianceK3[0,0]  
        
        # Gesamtanpassung Hier wird auch eine mehrfachbetrachtung des Offsets rausgerechnet (erster Term)
        summe = (paramsK1[3]+paramsK2[3]+paramsK3[3])*(-2/3)+(LorentzianOff(datav2[:,0],paramsK1[0],paramsK1[1],paramsK1[2],paramsK1[3])+LorentzianOff(datav2[:,0],paramsK2[0],paramsK2[1],paramsK2[2],paramsK2[3])+LorentzianOff(datav2[:,0],paramsK3[0],paramsK3[1],paramsK3[2],paramsK3[3]))
        ax1[Zeile,Spalte].plot(datav2[:,0],summe,color='red',linewidth=0.9)
        ax1[Zeile,Spalte].tick_params(axis='both', direction='in')

# Beschriftungen hinzufügen
ax1[0,0].set_ylabel('P0 (Bleistift)')
ax1[1,0].set_ylabel('P1')
ax1[2,0].set_ylabel('P2')
ax1[3,0].set_ylabel('P3')
ax1[3,0].set_xlabel('Ramanverschiebung (1/cm)')
ax1[3,1].set_xlabel('Ramanverschiebung (1/cm)')
ax1[0,0].set_title('parallele Polarisation')
ax1[0,1].set_title('senkrechte Polarisation')
ax1[1,0].set_yticks([0,5,10]) # Unterdrückt den Tick bei 15
ax1[2,0].text(500,25,'Intensität (arb. u.)',rotation=90)
    
# Berechnung der Depolarisationsraten rho = senkrecht/parallel
d = {'K1':[],'K1_err':[],'K2':[],'K2_err':[],'K3':[],'K3_err':[]}
depol = pd.DataFrame(data=d)

Messungen = [[0,5],[2,7],[3,8],[4,9]]
for i in range(4):
    parallel = Messungen[i][0]
    senkrecht = Messungen[i][1]
    
    depol.at[i,'K1'] =  graphitObject.at[senkrecht,'K1']/graphitObject.at[parallel,'K1']
    depol.at[i,'K1_err'] = depol.at[i,'K1']*np.sqrt((np.sqrt(graphitObject.at[senkrecht,'K1cov'])/graphitObject.at[senkrecht,'K1'])**2++(np.sqrt(graphitObject.at[parallel,'K1cov'])/graphitObject.at[parallel,'K1'])**2)
    depol.at[i,'K2'] =  graphitObject.at[senkrecht,'K2']/graphitObject.at[parallel,'K2']
    depol.at[i,'K2_err'] = depol.at[i,'K2']*np.sqrt((np.sqrt(graphitObject.at[senkrecht,'K2cov'])/graphitObject.at[senkrecht,'K2'])**2++(np.sqrt(graphitObject.at[parallel,'K2cov'])/graphitObject.at[parallel,'K2'])**2)
    depol.at[i,'K3'] =  graphitObject.at[senkrecht,'K3']/graphitObject.at[parallel,'K3']
    depol.at[i,'K3_err'] = depol.at[i,'K3']*np.sqrt((np.sqrt(graphitObject.at[senkrecht,'K3cov'])/graphitObject.at[senkrecht,'K3'])**2++(np.sqrt(graphitObject.at[parallel,'K3cov'])/graphitObject.at[parallel,'K3'])**2)

###-- Auswertung der kompletten unpolarisierten Graphit Spektren --##

# Betrachteter Bereich
min = 1000
max = 3000

# Fit-Bereiche für einzelne Peaks
minK1 = 1250
maxK1 = 1450
minK2 = 1400
maxK2 = 2250
minK3 = 2250
maxK3 = 2800

# Position und Breite der Peaks werden zur Zuordnung der zweiten Ordnung gesammelt
d = {'PosK1':[],'PosK1_err':[],'PosK2':[],'PosK2_err':[],'PosK3':[],'PosK3_err':[],'WidK1':[],'WidK1_err':[],'WidK2':[],'WidK2_err':[],'WidK3':[],'WidK3_err':[]}
zweiteOrdnung = pd.DataFrame(data=d)

fig1, ax1 = plt.subplots(4, 1, sharex ='col')
plt.subplots_adjust(wspace=0,hspace=0) # Kein Platz zwischen einzelnen Subplots
Messungen = [10,11,12,13]

for i in range(4):
    data_im = graphitObject.at[Messungen[i],'Daten']
    datav = np.transpose([data_im[:,0]+verschiebung,data_im[:,1]]) # Verschieben der Spektren um Wert aus Kalibriermessung
    datav2 = giveArea(datav,min,max) # Betrachtung eines Ausschnittes der Daten
    ax1[i].plot(datav2[:,0],datav2[:,1],'.',markersize=2,color='black') # Plotten der Messdaten

    # Anlegen der einzelnen Fits
    # K1
    ausschnittK1 = giveArea(datav2,minK1,maxK1) # Auschnitt zum Fitten
    paramsK1, params_covarianceK1 = optimize.curve_fit(LorentzianOff,ausschnittK1[:,0],ausschnittK1[:,1],p0=[0,1330,1,4]) # Fitten
    ax1[i].plot(datav2[:,0],LorentzianOff(datav2[:,0],paramsK1[0],paramsK1[1],paramsK1[2],paramsK1[3]),color='g',linewidth=1) # Plotten
    zweiteOrdnung.at[i,'PosK1'] = paramsK1[1]
    zweiteOrdnung.at[i,'PosK1_err'] = np.sqrt(abs(params_covarianceK1[1,1]))
    zweiteOrdnung.at[i,'WidK1'] = paramsK1[2]
    zweiteOrdnung.at[i,'WidK1_err'] = np.sqrt(abs(params_covarianceK1[2,2]))
    
    # K2
    ausschnittK2 = giveArea(datav2,minK2,maxK2)
    paramsK2, params_covarianceK2 = optimize.curve_fit(LorentzianOff,ausschnittK2[:,0],ausschnittK2[:,1],p0=[1,1575,1,4])
    ax1[i].plot(datav2[:,0],LorentzianOff(datav2[:,0],paramsK2[0],paramsK2[1],paramsK2[2],paramsK2[3]),color='g',linewidth=1)
    zweiteOrdnung.at[i,'PosK2'] = paramsK2[1]
    zweiteOrdnung.at[i,'PosK2_err'] = np.sqrt(abs(params_covarianceK2[1,1]))
    zweiteOrdnung.at[i,'WidK2'] = paramsK2[2]
    zweiteOrdnung.at[i,'WidK2_err'] = np.sqrt(abs(params_covarianceK2[2,2]))
    
    # K3
    ausschnittK3 = giveArea(datav2,minK3,maxK3)
    paramsK3, params_covarianceK3 = optimize.curve_fit(LorentzianOff,ausschnittK3[:,0],ausschnittK3[:,1],p0=[1,2650,1,4])
    ax1[i].plot(datav2[:,0],LorentzianOff(datav2[:,0],paramsK3[0],paramsK3[1],paramsK3[2],paramsK3[3]),color='g',linewidth=1)
    zweiteOrdnung.at[i,'PosK3'] = paramsK3[1]
    zweiteOrdnung.at[i,'PosK3_err'] = np.sqrt(abs(params_covarianceK3[1,1]))
    zweiteOrdnung.at[i,'WidK3'] = paramsK3[2]
    zweiteOrdnung.at[i,'WidK3_err'] = np.sqrt(abs(params_covarianceK3[2,2]))
    
    # Gesamtanpassung - Hier wird auch eine Mehrfachbetrachtung des Offsets rausgerechnet (erster Term)
    summe = (paramsK1[3]+paramsK2[3]+paramsK3[3])*(-2/3)+(LorentzianOff(datav2[:,0],paramsK1[0],paramsK1[1],paramsK1[2],paramsK1[3])+LorentzianOff(datav2[:,0],paramsK2[0],paramsK2[1],paramsK2[2],paramsK2[3])+LorentzianOff(datav2[:,0],paramsK3[0],paramsK3[1],paramsK3[2],paramsK3[3]))
    ax1[i].plot(datav2[:,0],summe,color='red',linewidth=0.9)
    ax1[i].tick_params(axis='both', direction='in')

# Beschriftung
ax1[0].set_ylabel('P0 (Bleistift)')
ax1[1].set_ylabel('P1')
ax1[2].set_ylabel('P2')
ax1[3].set_ylabel('P3')
ax1[3].set_xlabel('Ramanverschiebung (1/cm)')
ax1[2].text(625,100,'Intensität (arb. u.)',rotation=90)
    
##---- Vergleich von Spektren zweiter Ordnung ----#
daten = np.loadtxt('kalibrierung_allg.txt')
minS = 800
maxS = 1100
minG = 2500
maxG = 2800
fig1, ax = plt.subplots(2, 1)

Messungen = [10,11,12,13]

# 10
data_im = graphitObject.at[Messungen[0],'Daten']
datav = np.transpose([data_im[:,0]+verschiebung,data_im[:,1]]) # Verschieben der Spektren um Wert aus Kalibriermessung
datav2 = giveArea(datav,minG,maxG) # Betrachtung eines Ausschnittes der Daten
ax[1].plot(datav2[:,0],datav2[:,1],'.',color='black',markersize=2) # Plotten der Messdaten

# 11
data_im = graphitObject.at[Messungen[1],'Daten']
datav = np.transpose([data_im[:,0]+verschiebung,data_im[:,1]]) # Verschieben der Spektren um Wert aus Kalibriermessung
datav2 = giveArea(datav,minG,maxG) # Betrachtung eines Ausschnittes der Daten
ax[1].plot(datav2[:,0],1.5*datav2[:,1],'.',color='blue',markersize=3) # Plotten der Messdaten

# 12
data_im = graphitObject.at[Messungen[2],'Daten']
datav = np.transpose([data_im[:,0]+verschiebung,data_im[:,1]]) # Verschieben der Spektren um Wert aus Kalibriermessung
datav2 = giveArea(datav,minG,maxG) # Betrachtung eines Ausschnittes der Daten
ax[1].plot(datav2[:,0],10+datav2[:,1],'+',color='red') # Plotten der Messdaten

# 13
data_im = graphitObject.at[Messungen[3],'Daten']
datav = np.transpose([data_im[:,0]+verschiebung,data_im[:,1]]) # Verschieben der Spektren um Wert aus Kalibriermessung
datav2 = giveArea(datav,minG,maxG) # Betrachtung eines Ausschnittes der Daten
ax[1].plot(datav2[:,0],5+1.5*datav2[:,1],'x',color='green') # Plotten der Messdaten

# Silizium
ausschnittDarst = giveArea(daten,minS,maxS)
ausschnittFit = giveArea(daten,minS,maxS)
ax[0].plot(ausschnittDarst[:,0],ausschnittDarst[:,1],'.',color='black', markersize=2)

ax[0].set_ylabel('Silizium')
ax[1].set_ylabel('Graphit')
ax[1].set_xlabel('Ramanverschiebung (1/cm)')
ax[1].text(2460,85,'Intensität',rotation=90)
ax[1].legend(['P0','P1','P2','P3'],frameon=False)
ax[0].tick_params(axis='both', direction='in')
ax[1].tick_params(axis='both', direction='in')
ax[0].set_yticks([])
ax[1].set_yticks([])
plt.savefig('VergleichZweiterOrdnung.pdf')