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




filename = 'rotation_phase_diagram_gamma2_gamma1=1.25.txt'
beta11 = np.loadtxt(filename,delimiter='\t', skiprows=1, usecols = 0)
beta12 = np.loadtxt(filename,delimiter='\t', skiprows=1, usecols = 1)
beta21 = np.loadtxt(filename,delimiter='\t', skiprows=1, usecols = 2)
beta22 = np.loadtxt(filename,delimiter='\t', skiprows=1, usecols = 3)
ratev = np.loadtxt(filename,delimiter='\t', skiprows=1, usecols = 4)



print_sign = False

beta1 = 1/(0.5/beta21 + 0.5/beta11)#*np.sign(beta21)#*np.sign(beta11)
beta2 = 1/(0.5/beta22 + 0.5/beta12)#*np.sign(beta22)#*np.sign(beta12)



colormapp = plt.cm.RdBu_r

ttt = 0.05
bound = 1
ending = 150
starting = -150
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np
font = {'family': 'Arial', 'weight': 'normal', 'size': 10}    
NCURVES = n = len(ratev)
values = range(0, int(ending-starting), 1)
fig = plt.figure(figsize=(3.1, 3), dpi = 600)
ax = fig.add_subplot()
#ax.patch.set_facecolor('grey')
jet = cm = plt.get_cmap('jet') 
cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=colormapp)#coolwarm)#seismic)#jet)
print (scalarMap.get_clim())
lines = []
colorVal_rate = []
for i in range(0, n):
    if 0 <= int(ratev[i]*100 - starting) < int(ending-starting):
        colorVal_rate.append(scalarMap.to_rgba(values[int(ratev[i]*100 - starting)]))
    if int(ratev[i]*100 - starting) >= int(ending-starting):
        colorVal_rate.append(scalarMap.to_rgba(values[int(ending-starting)-1]))
    if int(ratev[i]*100 - starting) < 0:
        colorVal_rate.append(scalarMap.to_rgba(values[0]))
    
ax.scatter(beta1, beta2, color=colorVal_rate, s= 20, alpha = 0.9, marker = 'x', linewidth = 2.0)#,edgecolor = 'none')
ax.set_xticks([-8, -4, 0, 4, 8])
ax.set_yticks([-8, -4, 0, 4, 8])
ax.set_xticks([])
ax.set_yticks([])
ax=plt.gca();
ax.spines['bottom'].set_linewidth(ttt);
ax.spines['top'].set_linewidth(ttt);
ax.spines['left'].set_linewidth(ttt);
ax.spines['right'].set_linewidth(ttt);
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.legend(framealpha = 0.0)
plt.tick_params(labelsize = 10)
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Arial') for label in labels]
ax.legend(framealpha = 0, loc = 'lower right')
plt.tight_layout()
plt.show()

