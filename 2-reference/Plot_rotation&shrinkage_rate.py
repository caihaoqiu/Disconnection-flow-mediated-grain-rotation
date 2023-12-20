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


N = 300
Lx = 300
Ly = 300

for i in range(0, N+1):
    X.append(-Lx/2+Lx/N*i)
    Y.append(-Ly/2+Ly/N*i) 
x, y = np.meshgrid(X, Y)   


x1 = []
x2 = []
R = 100 
t_1 = []
area = []
w3_1 = []

All = 3000
Step =  int(All/50)+1

for iii in range(0, Step):
    t = 50*iii
    t_1.append(t)
    ratio = -1.0
    
    b1 = 1 
    b2 = 1.5 * b1
    
    h11 = 1*b2
    h21 = -2*b2
    h12 = -1*b1
    h22 = 2*b1

    
    beta11 = b1/h11
    beta21 = b1/h21
    beta12 = b2/h12
    beta22 = b2/h22

    M_rho1 = 0.01
    M_rho2 = 0.01
    path='./data/'
    x1_0 = np.loadtxt(path+'2mode-beta11_'+str(int(10*beta11))+'-beta12_'+str(int(10*beta12))+'-beta21_'+str(int(beta21*10))+'-beta22_'+str(int(beta22*10))+'-iniratio_'+str(ratio)+'-Mrho1_' + str(M_rho1)+'-Mrho2_' + str(M_rho2) +'-t_'+str(t) +'.txt',delimiter='\t', skiprows=1, usecols = 0)
    x2_0 = np.loadtxt(path+'2mode-beta11_'+str(int(10*beta11))+'-beta12_'+str(int(10*beta12))+'-beta21_'+str(int(beta21*10))+'-beta22_'+str(int(beta22*10))+'-iniratio_'+str(ratio)+'-Mrho1_' + str(M_rho1)+'-Mrho2_' + str(M_rho2) +'-t_'+str(t) +'.txt',delimiter='\t', skiprows=1, usecols = 1)
    area.append((abs(np.trapz(x2_0, x1_0)))/(np.pi*R**2))
    x1_0 = list(x1_0)
    x2_0 = list(x2_0)
    x1_0.append(x1_0[0])
    x2_0.append(x2_0[0])
    x1.append(x1_0)
    x2.append(x2_0)
    rho11 = np.loadtxt(path+'2mode-beta11_'+str(int(10*beta11))+'-beta12_'+str(int(10*beta12))+'-beta21_'+str(int(beta21*10))+'-beta22_'+str(int(beta22*10))+'-iniratio_'+str(ratio)+'-Mrho1_' + str(M_rho1)+'-Mrho2_' + str(M_rho2) +'-t_'+str(t) +'.txt',delimiter='\t', skiprows=1, usecols = 2)
    rho12 = np.loadtxt(path+'2mode-beta11_'+str(int(10*beta11))+'-beta12_'+str(int(10*beta12))+'-beta21_'+str(int(beta21*10))+'-beta22_'+str(int(beta22*10))+'-iniratio_'+str(ratio)+'-Mrho1_' + str(M_rho1)+'-Mrho2_' + str(M_rho2) +'-t_'+str(t) +'.txt',delimiter='\t', skiprows=1, usecols = 3)
    rho21 = np.loadtxt(path+'2mode-beta11_'+str(int(10*beta11))+'-beta12_'+str(int(10*beta12))+'-beta21_'+str(int(beta21*10))+'-beta22_'+str(int(beta22*10))+'-iniratio_'+str(ratio)+'-Mrho1_' + str(M_rho1)+'-Mrho2_' + str(M_rho2) +'-t_'+str(t) +'.txt',delimiter='\t', skiprows=1, usecols = 4)
    rho22 = np.loadtxt(path+'2mode-beta11_'+str(int(10*beta11))+'-beta12_'+str(int(10*beta12))+'-beta21_'+str(int(beta21*10))+'-beta22_'+str(int(beta22*10))+'-iniratio_'+str(ratio)+'-Mrho1_' + str(M_rho1)+'-Mrho2_' + str(M_rho2) +'-t_'+str(t) +'.txt',delimiter='\t', skiprows=1, usecols = 5)

    
    n = len(x1_0)-1
    
    
    a = 0.01
    G = 0.001
    
    kkk = 1
    


    mm = int(N/2)
    nn = int(N/2)
    J12_11 = 0
    J12_12 = 0
    J12_21 = 0
    J12_22 = 0        
     
    for j in range(0,n):                
        q1 = (X[mm]-x1[kkk-1][j])**2 + (Y[nn]-x2[kkk-1][j])**2 + a**2
        w12_1 = ((X[mm]-x1[kkk-1][j])/q1)
        w12_2 = ((Y[nn]-x2[kkk-1][j])/q1)           

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


        
        J12_11 += h11*rho11[j] * w12_1 *ll
        J12_12 += -h12*rho12[j] * w12_2 *ll
        J12_21 += h21*rho21[j] * w12_1 *ll
        J12_22 += -h22*rho22[j] * w12_2 *ll
        
                    
    w3 = 180/np.pi**2*(beta11*J12_11 + beta12*J12_12 + beta21*J12_21 + beta22*J12_22)
    w3_1.append(w3)
        

x_limit = 150 + 10
All = 4000
font = {'family': 'Arial', 'weight': 'normal', 'size': 10}
fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(x1[0], x2[0], c="red", linestyle = '--', linewidth = 2)
for iii in range (1,len(x1)):
#for iii in range (0,2):
        ax.plot(x1[iii], x2[iii], color = 'black', linewidth = 2)
#        ax.plot(rho1_1[iii], rho1_2[iii], linewidth = 2.0)
ax.set_xlim(xmin=-x_limit, xmax=x_limit)
ax.set_ylim(ymin=-x_limit, ymax=x_limit)
#ax.set_xlabel(x1, font)
#ax.set_ylabel(, font)
#ax.set_xticks([-100, -50, 0, 50, 100])
#ax.set_yticks([-100, -50, 0, 50, 100])
plt.tick_params(labelsize = 20)
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Arial') for label in labels]
plt.show()

font = {'family': 'Arial', 'weight': 'normal', 'size': 15}
fig, ax1 = plt.subplots(figsize=(4, 3), dpi = 600)
#ax1.scatter(t_1, sum_1, c="black")
ax1.plot(t_1, area, c='black', linewidth = 2.0)
ax1.set_xlim(xmin=0, xmax= All)
ax1.set_ylim(ymin=0, ymax=2)
#ax1.set_xticks([0,500, 1000, 1500, 2000, 2500])
#ax1.set_xticks([0,All/4, All/2, All*3/4, All])
#ax1.set_yticks([0,0.2, 0.4, 0.6, 0.8, 1.0])
ax1.set_xlabel(r't', font)
ax1.set_ylabel(r'[A(0)-A(t)]/A(0)', font)
plt.tick_params(labelsize = 16)
labels = ax1.get_xticklabels() + ax1.get_yticklabels()
[label.set_fontname('Arial') for label in labels]
#ax1.set_xlim(xmin=0, xmax=5100)
#ax1.set_ylim(ymin=0, ymax=100)
plt.show()


        
font = {'family': 'Arial', 'weight': 'normal', 'size': 16}
fig, ax = plt.subplots(figsize=(4,3), dpi =600)
ax.plot(t_1, w3_1, color = 'black')
ax.set_xlim(xmin=0, xmax=t)
ax.set_ylim(ymin=-4, ymax=4)
ax.set_xticks([0, 1000, 2000, 3000, 4000])
ax.set_xlabel(r'$t$', font)
ax.set_ylabel(r'rotation angle $\Delta \theta$', font)
plt.show()
