# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 17:13:17 2022

@author: USER
"""

import pandas as pd
import matplotlib.pyplot as plt
from numpy import log as ln
import math
import numpy as np

X=[]
Y=[]
F=[]

Kb=1.38e-23            #Boltzman Constant
T=288                  #absolute temperature
u=1.85e-5             #air dynamic viscocity
Uo=5e-2          #face velocity

dummy=0
f2=0

W=0.7e-3           #basic weight
rf=1220            #fiber material density
Z=11.8e-6          #fiber thickness
a=W/(rf*Z)      #alpha(fiber packing density)
#print(a)

df=208e-9       #mean fiber diameter  800 nm
Dp=40e-9           #particle diameter
x=0.1e-9           #pertambahan diameter partikel
lamda=67e-9         #min free path at ambient pressure
#Kn=2*lamda/Dp
#Knf=2*lamda/df
#print(Knf)

#c=Dp/df

#Ku=-(ln(a)/2)+a-(a**2)/4-0.75        #Kubawara hydrodunamic factor

#nr=0.6*((1-a)/Ku)*(1+(Knf/c))*(c**2/(1+c))

#penentuan nilai efisiensi D
#Cs=1+Kn*(1.207+0.44*math.exp(-0.78/Kn))     #cuningham slip
#D=(Kb*T*Cs)/(3*math.pi*u*Dp)                  #diffusion coefficient
#Pe=Uo*df/D                                  #Peclet number

#C1=1+0.388*Knf*((1-a)*Pe/Ku)**(1/3)
#C2=1/(1+1.6*((1-a)/Ku)**(1/3)*Pe**(-2/3)*C1)

#nd=1.6*((1-a)/Ku)**(1/3)*Pe**(-2/3)*C1*C2

#nf=nd+nr                        #correlation model

#cons=(4*a*nf*Z)/(math.pi*(1-a)*df)
#n=1-(math.exp(-(4*a*nf*Z)/(math.pi*(1-a)*df)))

while Dp<480e-9:
    Kn=2*lamda/Dp
    Knf=2*lamda/df
    c=Dp/df
    Ku=-(ln(a))/2+a-(a**2)/4-0.75        #Kubawara hydrodunamic factor

    nr=0.6*((1-a)/Ku)*(1+(Knf/c))*(c**2/(1+c))
    Cs=1+Kn*(1.207+0.44*math.exp(-0.78/Kn))     #cuningham slip
    D=Kb*T*Cs/(3*math.pi*u*Dp)                  #diffusion coefficient
    Pe=Uo*df/D                                  #Peclet number

    C1=1+0.388*Knf*((1-a)*Pe/Ku)**(1/3)
    C2=1/(1+1.6*((1-a)/Ku)**(1/3)*Pe**(-2/3)*C1)

    nd=1.6*((1-a)/Ku)**(1/3)*Pe**(-2/3)*C1*C2

    nf=nd+nr                        #correlation model

    n=1-(math.exp(-(4*a*nf*Z)/(math.pi*(1-a)*df)))
    f=(n-dummy)/x
    #simpan hasip Dp dan effisiensi
    X.append(Dp)
    Y.append(n)
    F.append(f)
    if (f2<0) and (f>0):
        print("nilai MPPS")
        print(Dp)
    
    #ukuran partikel ditambah dengan x
    Dp=Dp+x
    dummy=n         #dummy untuk menyimpan nilai eff di Dp sebelumnya
    f2=f            #untuk mencari nilai Dp

plt.plot(X,Y)
#print(F)




