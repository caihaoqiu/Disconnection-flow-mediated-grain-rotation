# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 13:04:58 2021

@author: qch1151
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
import math


X = []
Y = []


N = 100
Lx = 300
Ly = 300

for i in range(0, N+1):
    X.append(-Lx/2+Lx/N*i)
    Y.append(-Ly/2+Ly/N*i) 
x, y = np.meshgrid(X, Y)   


x1 = []
x2 = []


t = 0
sss = 29
T = 1.1
starting = 0
ending = 4000
numm = int(ending/50)+1

if sss == 29:
    ratio1 = 0.41379310344827586
    b1 = np.sqrt(29)/29
    b2 = b1*np.sqrt(2)
    b3 = b1 
    b4 = b2
    
    h11 = 6 * b3
    h21 = -8.5 * b3
    h12 = 6 * b4
    h22 = -8.5* b4
    h13 = - h11
    h23 = - h21
    h14 = - h12
    h24 = - h22
    
    beta110 = b1 / h11
    beta210 = b1 / h21
    beta120 = b2 / h12
    beta220 = b2 / h22
    beta130 = b3 / h13
    beta230 = b3 / h23
    beta140 = b4 / h14
    beta240 = b4 / h24   


R = 100
w3 = []
t_1 = []
area = []
x1 = []
x2 = []
Time = np.linspace(starting, ending, numm)
for tt in Time:
    t = int(tt)
    t_1.append(t)
    x1_0 = []
    x2_0 = []
    rho11 = []
    rho12 = []
    rho21 = []
    rho22 = []
    rho23 = []
    rho24 = []
    rho13 = []
    rho14 = []
    path = './data'
    x1_0 = np.loadtxt(str(path)+'/2mode-beta11_'+str(int(10*beta110))+'-beta12_'+str(int(10*beta120))+'-beta13_'+str(int(10*beta130))+'-beta14_'+str(int(10*beta140))+'-beta21_'+str(int(beta210*10))+'-beta22_'+str(int(beta220*10))+'-beta23_'+str(int(beta230*10))+'-beta24_'+str(int(beta240*10)) +'-iniratio_'+str(ratio1)  + '-T_' + str(T) +'-t_'+str(t) +'.txt',delimiter='\t', skiprows=1, usecols = 0)
    x2_0 = np.loadtxt(str(path)+'/2mode-beta11_'+str(int(10*beta110))+'-beta12_'+str(int(10*beta120))+'-beta13_'+str(int(10*beta130))+'-beta14_'+str(int(10*beta140))+'-beta21_'+str(int(beta210*10))+'-beta22_'+str(int(beta220*10))+'-beta23_'+str(int(beta230*10))+'-beta24_'+str(int(beta240*10)) +'-iniratio_'+str(ratio1)  + '-T_' + str(T) +'-t_'+str(t) +'.txt',delimiter='\t', skiprows=1, usecols = 1)
    area.append((abs(np.trapz(x2_0, x1_0)))/(np.pi*R**2))
    x1_0 = list(x1_0)
    x2_0 = list(x2_0)
    x1_0.append(x1_0[0])
    x2_0.append(x2_0[0])
    x1.append(x1_0)
    x2.append(x2_0)

    rho11 = np.loadtxt(str(path)+'/2mode-beta11_'+str(int(10*beta110))+'-beta12_'+str(int(10*beta120))+'-beta13_'+str(int(10*beta130))+'-beta14_'+str(int(10*beta140))+'-beta21_'+str(int(beta210*10))+'-beta22_'+str(int(beta220*10))+'-beta23_'+str(int(beta230*10))+'-beta24_'+str(int(beta240*10)) +'-iniratio_'+str(ratio1)  + '-T_' + str(T) +'-t_'+str(t) +'.txt',delimiter='\t', skiprows=1, usecols = 2)
    rho12 = np.loadtxt(str(path)+'/2mode-beta11_'+str(int(10*beta110))+'-beta12_'+str(int(10*beta120))+'-beta13_'+str(int(10*beta130))+'-beta14_'+str(int(10*beta140))+'-beta21_'+str(int(beta210*10))+'-beta22_'+str(int(beta220*10))+'-beta23_'+str(int(beta230*10))+'-beta24_'+str(int(beta240*10)) +'-iniratio_'+str(ratio1)  + '-T_' + str(T) +'-t_'+str(t) +'.txt',delimiter='\t', skiprows=1, usecols = 3)
    rho13 = np.loadtxt(str(path)+'/2mode-beta11_'+str(int(10*beta110))+'-beta12_'+str(int(10*beta120))+'-beta13_'+str(int(10*beta130))+'-beta14_'+str(int(10*beta140))+'-beta21_'+str(int(beta210*10))+'-beta22_'+str(int(beta220*10))+'-beta23_'+str(int(beta230*10))+'-beta24_'+str(int(beta240*10)) +'-iniratio_'+str(ratio1)  + '-T_' + str(T) +'-t_'+str(t) +'.txt',delimiter='\t', skiprows=1, usecols = 4)
    rho14 = np.loadtxt(str(path)+'/2mode-beta11_'+str(int(10*beta110))+'-beta12_'+str(int(10*beta120))+'-beta13_'+str(int(10*beta130))+'-beta14_'+str(int(10*beta140))+'-beta21_'+str(int(beta210*10))+'-beta22_'+str(int(beta220*10))+'-beta23_'+str(int(beta230*10))+'-beta24_'+str(int(beta240*10)) +'-iniratio_'+str(ratio1)  + '-T_' + str(T) +'-t_'+str(t) +'.txt',delimiter='\t', skiprows=1, usecols = 5)
    rho21 = np.loadtxt(str(path)+'/2mode-beta11_'+str(int(10*beta110))+'-beta12_'+str(int(10*beta120))+'-beta13_'+str(int(10*beta130))+'-beta14_'+str(int(10*beta140))+'-beta21_'+str(int(beta210*10))+'-beta22_'+str(int(beta220*10))+'-beta23_'+str(int(beta230*10))+'-beta24_'+str(int(beta240*10)) +'-iniratio_'+str(ratio1)  + '-T_' + str(T) +'-t_'+str(t) +'.txt',delimiter='\t', skiprows=1, usecols = 6)
    rho22 = np.loadtxt(str(path)+'/2mode-beta11_'+str(int(10*beta110))+'-beta12_'+str(int(10*beta120))+'-beta13_'+str(int(10*beta130))+'-beta14_'+str(int(10*beta140))+'-beta21_'+str(int(beta210*10))+'-beta22_'+str(int(beta220*10))+'-beta23_'+str(int(beta230*10))+'-beta24_'+str(int(beta240*10)) +'-iniratio_'+str(ratio1)  + '-T_' + str(T) +'-t_'+str(t) +'.txt',delimiter='\t', skiprows=1, usecols = 7)
    rho23 = np.loadtxt(str(path)+'/2mode-beta11_'+str(int(10*beta110))+'-beta12_'+str(int(10*beta120))+'-beta13_'+str(int(10*beta130))+'-beta14_'+str(int(10*beta140))+'-beta21_'+str(int(beta210*10))+'-beta22_'+str(int(beta220*10))+'-beta23_'+str(int(beta230*10))+'-beta24_'+str(int(beta240*10)) +'-iniratio_'+str(ratio1)  + '-T_' + str(T) +'-t_'+str(t) +'.txt',delimiter='\t', skiprows=1, usecols = 8)
    rho24 = np.loadtxt(str(path)+'/2mode-beta11_'+str(int(10*beta110))+'-beta12_'+str(int(10*beta120))+'-beta13_'+str(int(10*beta130))+'-beta14_'+str(int(10*beta140))+'-beta21_'+str(int(beta210*10))+'-beta22_'+str(int(beta220*10))+'-beta23_'+str(int(beta230*10))+'-beta24_'+str(int(beta240*10)) +'-iniratio_'+str(ratio1)  + '-T_' + str(T) +'-t_'+str(t) +'.txt',delimiter='\t', skiprows=1, usecols = 9)
    
    n = len(x1_0)-1
    
    
    
    e1 = [1,0]
    e2 = [0,1]
    phi_ref=[]
    vector = []
    ref = 4
    for i in range(0, ref):
        phi_ref.append(math.pi/ref*i)
    for i in range(0, ref):
        phi_ref.append(math.pi/ref*i)
    for i in range(0, ref):
        vector.append([math.cos(math.pi/(ref)*i), math.sin(math.pi/(ref)*i)])
    for i in range(0, ref):
        vector.append([math.cos(math.pi/(ref)*i), math.sin(math.pi/(ref)*i)])
    phi_ref.append(0)
    vector.append([1,0])
    
    
    dd = n / (2*ref)
    coeff = math.sqrt(2)/2
    a = 0.01
    
    kkk = 1
    coeff = math.sqrt(2)/2
    J11_12 = 0
    J12_12 = 0
    J13_12 = 0
    J14_12 = 0        

    J21_12 = 0
    J22_12 = 0
    J23_12 = 0
    J24_12 = 0        
    for j in range(0,n):    
        i = int(N/2)
        k = int(N/2)
        q1 = (X[i]-x1[kkk-1][j])**2 + (Y[k]-x2[kkk-1][j])**2 + a**2

        w12_1 = ((X[i]-x1[kkk-1][j])/q1)
        w12_3 = -0.5*((Y[k]-x2[kkk-1][j])/q1) 
        w12_2 = coeff*(w12_1 + w12_3)
        w12_4 = coeff*(-w12_1 + w12_3) 
        if j == 0 : 
            devj1_1 = (x1[kkk-1][j+1]-x1[kkk-1][n-1])/2
            devj1_2 = (x2[kkk-1][j+1]-x2[kkk-1][n-1])/2  
        if j == n-1 : 
            devj1_1 = (x1[kkk-1][0]-x1[kkk-1][j-1])/2
            devj1_2 = (x2[kkk-1][0]-x2[kkk-1][j-1])/2        
        if (j != 0 and j != n-1): 
            devj1_1 = (x1[kkk-1][j+1]-x1[kkk-1][j-1])/2
            devj1_2 = (x2[kkk-1][j+1]-x2[kkk-1][j-1])/2
        ll = np.sqrt(devj1_1**2+devj1_2**2)            
        
        J11_12 += h11*rho11[j]*ll * w12_1
        J12_12 += h12*rho12[j]*ll * w12_2
        J21_12 += h21*rho21[j]*ll * w12_1
        J22_12 += h22*rho22[j]*ll * w12_2
        J13_12 += h13*rho13[j]*ll * w12_3
        J23_12 += h23*rho23[j]*ll * w12_3
        J14_12 += h14*rho14[j]*ll * w12_4
        J24_12 += h24*rho24[j]*ll * w12_4
        phi_ref1 = math.pi/ref * 0
        phi_ref2 = math.pi/ref * 1
        beta11 = 1/(math.sin(phi_ref2-phi_ref1)) * beta110
        beta12 = 1/(math.sin(phi_ref2-phi_ref1)) * beta120
        beta13 = 1/(math.sin(phi_ref2-phi_ref1)) * beta130
        beta14 = 1/(math.sin(phi_ref2-phi_ref1)) * beta140 
        beta21 = 1/(math.sin(phi_ref2-phi_ref1)) * beta210
        beta22 = 1/(math.sin(phi_ref2-phi_ref1)) * beta220
        beta23 = 1/(math.sin(phi_ref2-phi_ref1)) * beta230
        beta24 = 1/(math.sin(phi_ref2-phi_ref1)) * beta240 

        w3_l = 180/np.pi**2*(beta11*J11_12 + beta12*J12_12 + beta21*J21_12 + beta22*J22_12 + beta13*J13_12 + beta14*J14_12 + beta23*J23_12 + beta24*J24_12)
    
    if sss == 5:
        w3_l = w3_l + 36.9
    if sss == 17:
        w3_l = w3_l + 28
    if sss == 29:
        w3_l = w3_l + 44
    if sss == 25:
        w3_l = w3_l + 16


    w3.append(w3_l)

t_dimension = 1.539*1e-4

#MD_data_rotation
T_29 = np.loadtxt('sigma29-44-MD.txt', delimiter=',',usecols=0)
W_29 = np.loadtxt('sigma29-44-MD.txt', delimiter=',',usecols=1)

font = {'family': 'Arial', 'weight': 'normal', 'size': 15}
fig, ax = plt.subplots(figsize=(4,3), dpi =600)
ax.plot(T_29, W_29-W_29[0], color = 'blue', linewidth = 4.0, linestyle = ":", label = 'MD_data')
ax.plot(np.array(t_1)*t_dimension, w3-w3[0], color = 'black', label = 'Continuum_data')
ax.set_xlim(xmin=0, xmax=0.6)
ax.set_ylim(ymin=-2, ymax=5)
ax.set_xlabel(r'$t$ (ns)', font)
ax.set_ylabel(r'$\Delta \theta~(^\circ)$', font)
plt.legend()
plt.show() 


#MD_data_shrinkage
TT_29 = np.loadtxt('sigma29-a.txt', delimiter=',',usecols=0)
A_29 = np.loadtxt('sigma29-a.txt', delimiter=',',usecols=1)
A_29 = A_29/A_29[0]


font = {'family': 'Arial', 'weight': 'normal', 'size': 15}
fig, ax1 = plt.subplots(figsize=(4, 3), dpi = 600)
ax1.plot(TT_29, A_29, c='blue', linewidth = 4.0, linestyle = ":", label = 'MD_data')
ax1.plot(np.array(t_1)*t_dimension, area, c='black', linewidth = 2.0, label = 'Continuum_data')
ax1.set_xlim(xmin=0, xmax= 0.6)
ax1.set_ylim(ymin=0.5, ymax=1)
ax1.set_xlabel(r'$t$ (ns)', font)
ax1.set_ylabel(r'$[A(0)-A(t)]/A(0)$', font)
#plt.tick_params(labelsize = 16)
#labels = ax1.get_xticklabels() + ax1.get_yticklabels()
#[label.set_fontname('Arial') for label in labels]
plt.legend()
plt.show()

x_limit = 150 + 10
All = 4000
font = {'family': 'Arial', 'weight': 'normal', 'size': 10}
fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(x1[0], x2[0], c="red", linestyle = '--', linewidth = 2)
for iii in range (1,len(x1)):
    if iii % 20 == 0:
        ax.plot(x1[iii], x2[iii], color = 'black', linewidth = 2)
ax.set_xlim(xmin=-x_limit, xmax=x_limit)
ax.set_ylim(ymin=-x_limit, ymax=x_limit)
plt.tick_params(labelsize = 20)
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Arial') for label in labels]
plt.show()