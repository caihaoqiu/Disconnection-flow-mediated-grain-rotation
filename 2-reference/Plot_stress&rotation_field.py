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


N = 150
Lx = 300
Ly = 300

for i in range(0, N+1):
    X.append(-Lx/2+Lx/N*i)
    Y.append(-Ly/2+Ly/N*i) 
x, y = np.meshgrid(X, Y)   


x1 = []
x2 = []

#t = 0

t = 0
ratio = -1.0

b1 = 1
b2 = 1.5

h11 = b2
h21 = -2*b2
h12 = -b1
h22 = 2*b1


beta11 = b1/h11
beta21 = b1/h21
beta12 = b2/h12
beta22 = b2/h22

M_rho1 = 0.01
M_rho2 = 0.01


va = 3
va_s = 1.0*1e-6
num_iso = 100
    
Time = np.linspace(2000, 2000, 1)
for tt in Time:
    t = int(tt)
    x1_0 = []
    x2_0 = []
    x1 = []
    x2 = []
    rho11 = []
    rho12 = []
    rho21 = []
    rho22 = []
    path='./data/'
    x1_0 = np.loadtxt(path+'2mode-beta11_'+str(int(10*beta11))+'-beta12_'+str(int(10*beta12))+'-beta21_'+str(int(beta21*10))+'-beta22_'+str(int(beta22*10))+'-iniratio_'+str(ratio)+'-Mrho1_' + str(M_rho1)+'-Mrho2_' + str(M_rho2) +'-t_'+str(t) +'.txt',delimiter='\t', skiprows=1, usecols = 0)
    x2_0 = np.loadtxt(path+'2mode-beta11_'+str(int(10*beta11))+'-beta12_'+str(int(10*beta12))+'-beta21_'+str(int(beta21*10))+'-beta22_'+str(int(beta22*10))+'-iniratio_'+str(ratio)+'-Mrho1_' + str(M_rho1)+'-Mrho2_' + str(M_rho2) +'-t_'+str(t) +'.txt',delimiter='\t', skiprows=1, usecols = 1)
    x1.append(x1_0)
    x2.append(x2_0)
    rho11 = np.loadtxt(path+'2mode-beta11_'+str(int(10*beta11))+'-beta12_'+str(int(10*beta12))+'-beta21_'+str(int(beta21*10))+'-beta22_'+str(int(beta22*10))+'-iniratio_'+str(ratio)+'-Mrho1_' + str(M_rho1)+'-Mrho2_' + str(M_rho2) +'-t_'+str(t) +'.txt',delimiter='\t', skiprows=1, usecols = 2)
    rho12 = np.loadtxt(path+'2mode-beta11_'+str(int(10*beta11))+'-beta12_'+str(int(10*beta12))+'-beta21_'+str(int(beta21*10))+'-beta22_'+str(int(beta22*10))+'-iniratio_'+str(ratio)+'-Mrho1_' + str(M_rho1)+'-Mrho2_' + str(M_rho2) +'-t_'+str(t) +'.txt',delimiter='\t', skiprows=1, usecols = 3)
    rho21 = np.loadtxt(path+'2mode-beta11_'+str(int(10*beta11))+'-beta12_'+str(int(10*beta12))+'-beta21_'+str(int(beta21*10))+'-beta22_'+str(int(beta22*10))+'-iniratio_'+str(ratio)+'-Mrho1_' + str(M_rho1)+'-Mrho2_' + str(M_rho2) +'-t_'+str(t) +'.txt',delimiter='\t', skiprows=1, usecols = 4)
    rho22 = np.loadtxt(path+'2mode-beta11_'+str(int(10*beta11))+'-beta12_'+str(int(10*beta12))+'-beta21_'+str(int(beta21*10))+'-beta22_'+str(int(beta22*10))+'-iniratio_'+str(ratio)+'-Mrho1_' + str(M_rho1)+'-Mrho2_' + str(M_rho2) +'-t_'+str(t) +'.txt',delimiter='\t', skiprows=1, usecols = 5)
    
    n = len(x1_0)
    
    
    coeff = math.sqrt(2)/2
    a = 0.01
    G = 1e-5
    
    kkk = 1
    coeff = math.sqrt(2)/2
    
    I1_12 = I2_12 = I3_12 = I4_12 = 0        
    I1_11 = I2_11 = I3_11 = I4_11 = 0 
    I1_22 = I2_22 = I3_22 = I4_22 = 0 
    
    sigma11 = []
    sigma12 = []
    sigma22 = []
    w3 = []
    rotation = []
    
    

    for i in range (0, N+1):
        sigma11_l = []
        sigma12_l = []
        sigma22_l = []
        w3_l = []
        for k in range(0, N+1):
            I12_11 = 0
            I12_12 = 0
            I12_21 = 0
            I12_22 = 0
            J12_11 = 0
            J12_12 = 0
            J12_21 = 0
            J12_22 = 0        
            for j in range(0,n):                
                q1 = (X[i]-x1[kkk-1][j])**2 + (Y[k]-x2[kkk-1][j])**2 + a**2
                tao12_1 = ((X[i]-x1[kkk-1][j])/q1)*(1-2*(Y[k]-x2[kkk-1][j])**2/q1)    
                tao12_2 = ((Y[k]-x2[kkk-1][j])/q1)*(1-2*(Y[k]-x2[kkk-1][j])**2/q1)
                w12_1 = ((X[i]-x1[kkk-1][j])/q1)
                w12_2 = ((Y[k]-x2[kkk-1][j])/q1)           
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
                
                I12_11 += h11*rho11[j]*ll * tao12_1
                I12_12 += h12*rho12[j]*ll * tao12_2
                I12_21 += h21*rho21[j]*ll * tao12_1
                I12_22 += h22*rho22[j]*ll * tao12_2
                
                J12_11 += h11*rho11[j] * w12_1 *ll
                J12_12 += -h12*rho12[j] * w12_2 *ll
                J12_21 += h21*rho21[j] * w12_1 *ll
                J12_22 += -h22*rho22[j] * w12_2 *ll
                           
            ssss = G*(beta11*I12_11 + beta12*I12_12 + beta21*I12_21 + beta22*I12_22)
            rrrr = 180/np.pi**2*(beta11*J12_11 + beta12*J12_12 + beta21*J12_21 + beta22*J12_22)
            
            if ssss > 1*va_s: ssss = va_s
            if ssss < -1*va_s: ssss = -va_s
            if rrrr > 1*va: rrrr = va
            if rrrr < -1*va: rrrr = -va
            
            sigma12_l.append(np.round(ssss, 6))
            w3_l.append(np.round(rrrr, 3))
        sigma12.append(sigma12_l)
        w3.append(w3_l)
    '''
    if t == 0:
        w3_0 = w3
    
    rotation = np.array(w3) - np.array(w3_0)
    '''
#    print(np.size(rotation), np.size(w3))          
#    print(rotation)
        
    x1_0 = list(x1_0)
    x1_0.append(x1_0[0])
    x2_0 = list(x2_0)
    x2_0.append(x2_0[0])
    
    #print(x1_0)
    
    colorm = plt.cm.RdBu.reversed()#bwr#seismic#jet
    ttt = 0.05

    font = {'family': 'Arial', 'weight': 'normal', 'size': 16}
    fig, ax = plt.subplots(figsize=(3,3), dpi =600)
    ro = ax.contourf(y, x, sigma12, num_iso, vmin=-va_s, vmax=va_s, cmap=colorm)#seismic)#jet)
    ax.set_xlim(xmin=-150, xmax=150)
    ax.set_ylim(ymin=-150, ymax=150)
    ax.set_xticks([])
    ax.set_yticks([])
    ax=plt.gca();
    ax.spines['bottom'].set_linewidth(ttt);
    ax.spines['top'].set_linewidth(ttt);
    ax.spines['left'].set_linewidth(ttt);
    ax.spines['right'].set_linewidth(ttt);
    plt.tight_layout()
    plt.show()
    

    font = {'family': 'Arial', 'weight': 'normal', 'size': 16}
    fig, ax = plt.subplots(figsize=(3,3), dpi =600)
    plt.contourf(y, x, w3, num_iso, vmin=-va, vmax=va, cmap=colorm)#seismic)#jet)
    ax.set_xlim(xmin=-150, xmax=150)
    ax.set_ylim(ymin=-150, ymax=150)
    ax.set_xticks([])
    ax.set_yticks([])
    ax=plt.gca();
    ax.spines['bottom'].set_linewidth(ttt);
    ax.spines['top'].set_linewidth(ttt);
    ax.spines['left'].set_linewidth(ttt);
    ax.spines['right'].set_linewidth(ttt);
    plt.tight_layout()
    plt.show()