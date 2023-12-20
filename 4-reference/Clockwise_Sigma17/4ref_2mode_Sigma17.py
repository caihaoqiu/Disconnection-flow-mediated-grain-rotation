
import numpy as np
import matplotlib.pyplot as plt
import math
import os


def mo(x1, x2):
    length = math.sqrt(x1**2+x2**2)
    return length

def area(x, y):
    sum1 = 0
    for i in range(1, len(x)):
        slope = (x[i]-x[i-1])/(y[i] - y[i-1])
        delta_x = x[i] - x[i - 1]
        delta_y =  y[i - 1] + (y[i-1] + delta_x * slope)
        sum1 = sum1 + (delta_y/2 * delta_x)
    return sum1

def dev(step, N, Nmax): 
    if i == 0 : 
        qta = mo(x1[step-1][N+1]-x1[step-1][N],x2[step-1][N+1]-x2[step-1][N])-mo(x1[step-1][Nmax-1]-x1[step-1][N],x2[step-1][Nmax-1]-x2[step-1][N])
        dev1_1 = (x1[step-1][N+1]-x1[step-1][Nmax-1])/2
        dev1_2 = (x2[step-1][N+1]-x2[step-1][Nmax-1])/2
        dev2_1 = x1[step-1][N+1]+x1[step-1][Nmax-1]-2*x1[step-1][N]
        dev2_2 = x2[step-1][N+1]+x2[step-1][Nmax-1]-2*x2[step-1][N]        
    if i == Nmax-1 : 
        qta = mo(x1[step-1][0]-x1[step-1][N],x2[step-1][0]-x2[step-1][N])-mo(x1[step-1][N-1]-x1[step-1][N],x2[step-1][N-1]-x2[step-1][N])      
        dev1_1 = (x1[step-1][0]-x1[step-1][N-1])/2
        dev1_2 = (x2[step-1][0]-x2[step-1][N-1])/2
        dev2_1 = x1[step-1][0]+x1[step-1][N-1]-2*x1[step-1][N]
        dev2_2 = x2[step-1][0]+x2[step-1][N-1]-2*x2[step-1][N]    
    if (i != 0 and i != Nmax-1): 
        qta = mo(x1[step-1][N+1]-x1[step-1][N],x2[step-1][N+1]-x2[step-1][N])-mo(x1[step-1][N-1]-x1[step-1][N],x2[step-1][N-1]-x2[step-1][N])
        dev1_1 = (x1[step-1][N+1]-x1[step-1][N-1])/2
        dev1_2 = (x2[step-1][N+1]-x2[step-1][N-1])/2
        dev2_1 = x1[step-1][N+1]+x1[step-1][N-1]-2*x1[step-1][N]
        dev2_2 = x2[step-1][N+1]+x2[step-1][N-1]-2*x2[step-1][N]
    return(qta, dev1_1, dev1_2, dev2_1, dev2_2)
def devj(step, N, Nmax): 
    if j == 0 : 
        qta = mo(x1[step-1][N+1]-x1[step-1][N],x2[step-1][N+1]-x2[step-1][N])-mo(x1[step-1][Nmax-1]-x1[step-1][N],x2[step-1][Nmax-1]-x2[step-1][N])
        dev1_1 = (x1[step-1][N+1]-x1[step-1][Nmax-1])/2
        dev1_2 = (x2[step-1][N+1]-x2[step-1][Nmax-1])/2
        dev2_1 = x1[step-1][N+1]+x1[step-1][Nmax-1]-2*x1[step-1][N]
        dev2_2 = x2[step-1][N+1]+x2[step-1][Nmax-1]-2*x2[step-1][N]        
    if j == Nmax-1 : 
        qta = mo(x1[step-1][0]-x1[step-1][N],x2[step-1][0]-x2[step-1][N])-mo(x1[step-1][N-1]-x1[step-1][N],x2[step-1][N-1]-x2[step-1][N])      
        dev1_1 = (x1[step-1][0]-x1[step-1][N-1])/2
        dev1_2 = (x2[step-1][0]-x2[step-1][N-1])/2
        dev2_1 = x1[step-1][0]+x1[step-1][N-1]-2*x1[step-1][N]
        dev2_2 = x2[step-1][0]+x2[step-1][N-1]-2*x2[step-1][N]    
    if (j != 0 and j != Nmax-1): 
        qta = mo(x1[step-1][N+1]-x1[step-1][N],x2[step-1][N+1]-x2[step-1][N])-mo(x1[step-1][N-1]-x1[step-1][N],x2[step-1][N-1]-x2[step-1][N])
        dev1_1 = (x1[step-1][N+1]-x1[step-1][N-1])/2
        dev1_2 = (x2[step-1][N+1]-x2[step-1][N-1])/2
        dev2_1 = x1[step-1][N+1]+x1[step-1][N-1]-2*x1[step-1][N]
        dev2_2 = x2[step-1][N+1]+x2[step-1][N-1]-2*x2[step-1][N]
    return(qta, dev1_1, dev1_2, dev2_1, dev2_2)
def stablize_l(ll):
    return 6*np.tanh(ll/6)
def stablize(ll):
    ee = 0.01
    return ee*np.tanh(ll/ee)
def Rn(a):
    eta = 1.3
    return eta/5*(np.log(2)+np.log(1+np.cosh(5*a/eta)))            




e1 = [1,0]
e2 = [0,1]

t0 = math.pi/180
n = 60
R0 = 100

delt = 0.001
K_sp = 0.01
tao_ext = 0
G = 1e-2
a = 0.01
Ml = -1000
ep = 0.1


All = 4000
Total_step = int(All/delt)


gamma1 = gamma3 = 0.91
gamma2 = gamma4 = 0.85

b1 = np.sqrt(17)/17
b2 = b1*np.sqrt(2)
b3 = b1 
b4 = b2

h11 = -2 * b3
h21 = 6.5 * b3
h12 = -2 * b4
h22 = 6.5* b4
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



#Dimensionless temperature
T = 1.38e-3 * 800

#Dimensionless activation energy
E11 = 1.4*1.0
E21 = 1.4*1.0
E12 = 1.4*0.25
E22 = 1.4*1.2
ratio_1 = np.exp((E11 - E21) /T)
ratio_2 = np.exp((E12 - E22) /T)
ratio_12 = np.exp((E11 - E12) /T)
print(ratio_1, ratio_2, ratio_12)
M_11 = 1
M_21 = M_11 * ratio_1
M_12 = M_11 * ratio_12
M_22 = M_12 * ratio_2
normalize = max(M_11, M_21, M_12, M_22)
M_11 = M_11 / normalize
M_21 = M_21 / normalize
M_12 = M_12 / normalize
M_22 = M_22 / normalize
rho_x_ratio = 0.1
M_rho11 = M_11 * rho_x_ratio
M_rho21 = M_21 * rho_x_ratio
M_rho12 = M_12 * rho_x_ratio
M_rho22 = M_22 * rho_x_ratio



Ra = 5*abs(beta110*h11)

theta = []
x1 = []
x2 = []
x1_0 = []
x2_0 = []
t = [0]
t_1 = [0]
sum_1 = []
x1_zero = []
for ii in range(0, n):
    theta.append(360/n*ii*t0)
    x1_0.append(R0*math.cos(theta[ii]))
    x2_0.append(R0*math.sin(theta[ii]))
x1_0.append(R0)
x2_0.append(0)
x1.append(x1_0)
x2.append(x2_0)
x1_zero.append(x1_0[0])
sum_1.append(abs(np.trapz(x2_0, x1_0)/R0**2))
phi_ref=[]
vectorx = [[1,0]]
vectory = []
vector = []
ref = 4
    
for i in range(0, ref):
    vector.append([math.cos(math.pi/(ref)*i), math.sin(math.pi/(ref)*i)])
for i in range(0, ref):
    vector.append([math.cos(math.pi/(ref)*i), math.sin(math.pi/(ref)*i)])

vector.append([1,0])

dd = n / (2*ref)
coeff = math.sqrt(2)/2
coeff1 = -math.sqrt(2)/2

l1_0 = []
l2_0 = []
rho1_1_0 = []
rho1_2_0 = []
rho2_1_0 = []
rho2_2_0 = []
rho1_3_0 = []
rho1_4_0 = []
rho2_3_0 = []
rho2_4_0 = []
length_0 = []
lamda_0 = []
ratio1 = 4/17
ratio2 = 4/17
ratio3 = 4/17
ratio4 = 4/17
sca = 1e-5
G_r = 1e-5
LL_sca = 2e-6
for i in range(0, n):
    dev1_1 = dev(1, i, n)[1]
    dev1_2 = dev(1, i, n)[2]
    dev2_1 = dev(1, i, n)[3]
    dev2_2 = dev(1, i, n)[4]    
    ll1 = dev1_1/np.sqrt(dev1_1**2+dev1_2**2)
    ll2 = dev1_2/np.sqrt(dev1_1**2+dev1_2**2)
    if ll1 != 0 and ll2/ll1 >= 0 and ll2/ll1 < 1:
        rho11 = ratio1*np.sqrt(2)*ll2 / h11
        rho21 = (1-ratio1)*np.sqrt(2)*ll2 / h21
        rho12 = -ratio2*(ll1 - ll2)/h12
        rho22 = -(1-ratio2)*(ll1 - ll2)/h22
        rho13 = rho23 = 0
        rho14 = rho24 = 0
    if ll1 != 0 and ll2/ll1 >= 1:
        rho12 = ratio2*(-ll1 + ll2)/h12
        rho22 = (1-ratio2)*(-ll1 + ll2)/h22
        rho13 = -ratio3*np.sqrt(2)*ll1 / h13
        rho23 = -(1-ratio3)*np.sqrt(2)*ll1 / h23
        rho11 = rho21 = 0
        rho14 = rho24 = 0
    if ll1 == 0 or (ll1 != 0 and ll2/ll1 < -1):       
        rho13 = -ratio3*np.sqrt(2)*ll1 / h13
        rho23 = -(1-ratio3)*np.sqrt(2)*ll1 / h23
        rho14 = -ratio4*(ll1 + ll2)/h14
        rho24 = -(1-ratio4)*(ll1 + ll2)/h24
        rho11 = rho21 = 0
        rho12 = rho22 = 0
    if ll1 != 0 and ll2/ll1 >= -1 and ll2/ll1 < 0:                
        rho11 = -ratio1*np.sqrt(2)*ll2 / h11
        rho21 = -(1-ratio1)*np.sqrt(2)*ll2 / h21
        rho14 = ratio4*(ll1 + ll2)/h14
        rho24 = (1-ratio4)*(ll1 + ll2)/h24
        rho13 = rho23 = 0
        rho12 = rho22 = 0
    l1_0.append(ll1)
    l2_0.append(ll2)       
    rho1_1_0.append(rho11)        #\rho_1^{(1)}
    rho1_2_0.append(rho12)        #\rho_1^{(2)}    
    rho2_1_0.append(rho21)        #\rho_2^{(1)}
    rho2_2_0.append(rho22)        #\rho_2^{(2)}    
    rho1_3_0.append(rho13)        #\rho_1^{(1)}
    rho1_4_0.append(rho14)        #\rho_1^{(2)}    
    rho2_3_0.append(rho23)        #\rho_2^{(1)}
    rho2_4_0.append(rho24)        #\rho_2^{(2)}    
    length_00 = np.sqrt(dev1_1**2 + dev1_2**2)
    length_0.append(length_00)
    lamda_0.append(0)

length = [length_0]
l1 = [l1_0]
l2 = [l2_0]

rho1_1 = [rho1_1_0]
rho1_2 = [rho1_2_0]
rho2_1 = [rho2_1_0]
rho2_2 = [rho2_2_0]
rho1_3 = [rho1_3_0]
rho1_4 = [rho1_4_0]
rho2_3 = [rho2_3_0]
rho2_4 = [rho2_4_0]
lamda = [lamda_0]



for kkk in range(1,Total_step+1):
    t.append(delt*kkk)
    x1_new = []
    x2_new = []
    tao1 = []  
    lamda_new = []
    rho1_1_new = []
    rho1_2_new = []
    rho2_1_new = []
    rho2_2_new = []
    rho1_3_new = []
    rho1_4_new = []
    rho2_3_new = []
    rho2_4_new = []
    length_new = []

    close = 0
    close_num = []
    for i in range(0, n):
        I11_12 = 0
        I12_12 = 0
        I13_12 = 0
        I14_12 = 0        
        I11_11 = 0
        I12_11 = 0
        I13_11 = 0
        I14_11 = 0        
        I11_22 = 0
        I12_22 = 0
        I13_22 = 0
        I14_22 = 0 

        I21_12 = 0
        I22_12 = 0
        I23_12 = 0
        I24_12 = 0        
        I21_11 = 0
        I22_11 = 0
        I23_11 = 0
        I24_11 = 0        
        I21_22 = 0
        I22_22 = 0
        I23_22 = 0
        I24_22 = 0 
        
        f_R_1_1 = 0
        f_R_2_1 = 0
        f_R_1_2 = 0
        f_R_2_2 = 0
        f_R_1_3 = 0
        f_R_2_3 = 0
        f_R_1_4 = 0
        f_R_2_4 = 0     
        
        f_th_1_1 = 0
        f_th_2_1 = 0
        f_th_1_2 = 0
        f_th_2_2 = 0
        f_th_1_3 = 0
        f_th_2_3 = 0
        f_th_1_4 = 0
        f_th_2_4 = 0
        
        qta = dev(kkk, i, n)[0]
        dev1_1 = dev(kkk, i, n)[1]
        dev1_2 = dev(kkk, i, n)[2]
        dev2_1 = dev(kkk, i, n)[3]
        dev2_2 = dev(kkk, i, n)[4]   

        l1_car = dev1_1/math.sqrt(dev1_1**2+dev1_2**2)
        l2_car = dev1_2/math.sqrt(dev1_1**2+dev1_2**2)     
        l = [l1_car, l2_car]
        l_ = [-l1_car, -l2_car]
        
        length_kkk = np.sqrt(dev1_1**2 + dev1_2**2)
        if l1_car != 0 and l2_car/l1_car >= 0 and l2_car/l1_car < 1:
            l1_old = - (h12 * rho1_2[kkk-1][i] + h22 * rho2_2[kkk-1][i]) * np.cos(0) + (h11 * rho1_1[kkk-1][i] + h21 * rho2_1[kkk-1][i]) * np.cos(np.pi/4)
            l2_old = - (h12 * rho1_2[kkk-1][i] + h22 * rho2_2[kkk-1][i]) * np.sin(0) + (h11 * rho1_1[kkk-1][i] + h21 * rho2_1[kkk-1][i]) * np.sin(np.pi/4)
        if l1_car != 0 and l2_car/l1_car >= 1:
            l1_old = - (h13 * rho1_3[kkk-1][i] + h23 * rho2_3[kkk-1][i]) * np.cos(np.pi/4) + (h12 * rho1_2[kkk-1][i] + h22 * rho2_2[kkk-1][i]) * np.cos(np.pi/2)
            l2_old = - (h13 * rho1_3[kkk-1][i] + h23 * rho2_3[kkk-1][i]) * np.sin(np.pi/4) + (h12 * rho1_2[kkk-1][i] + h22 * rho2_2[kkk-1][i]) * np.sin(np.pi/2)
        if l1_car == 0 or (l1_car != 0 and l2_car/l1_car < -1):       
            l1_old = - (h14 * rho1_4[kkk-1][i] + h24 * rho2_4[kkk-1][i]) * np.cos(np.pi/2) + (h13 * rho1_3[kkk-1][i] + h23 * rho2_3[kkk-1][i]) * np.cos(3*np.pi/4)
            l2_old = - (h14 * rho1_4[kkk-1][i] + h24 * rho2_4[kkk-1][i]) * np.sin(np.pi/2) + (h13 * rho1_3[kkk-1][i] + h23 * rho2_3[kkk-1][i]) * np.sin(3*np.pi/4)
        if l1_car != 0 and l2_car/l1_car >= -1 and l2_car/l1_car < 0:                
            l1_old =  -(h11 * rho1_1[kkk-1][i] + h21 * rho2_1[kkk-1][i]) * np.cos(3*np.pi/4) + (h14 * rho1_4[kkk-1][i] + h24 * rho2_4[kkk-1][i]) * np.cos(np.pi*0)
            l2_old =  -(h11 * rho1_1[kkk-1][i] + h21 * rho2_1[kkk-1][i]) * np.sin(3*np.pi/4) + (h14 * rho1_4[kkk-1][i] + h24 * rho2_4[kkk-1][i]) * np.sin(np.pi*0)
        lamda_new_value = lamda[kkk-1][i] + Ml*delt/length_kkk * (l1_old**2 + l2_old**2 - 1)
        lamda_new.append(lamda_new_value)




        for j in range(0,n):
            dx1 = x1[kkk-1][i]-x1[kkk-1][j]
            dx2 = x2[kkk-1][i]-x2[kkk-1][j]
            RR = np.sqrt(dx1**2 + dx2**2)
            if i != j:
                costheta1 = dx2/RR
                costheta2 = (-coeff * dx1 + coeff * dx2)/RR
                costheta3 = -dx1/RR
                costheta4 = (-coeff * dx1 - coeff * dx2)/RR
            dev1_1jjj = devj(kkk, j, n)[1]
            dev1_2jjj = devj(kkk, j, n)[2]
            length_jjj = np.sqrt(dev1_1jjj**2 + dev1_2jjj**2)


            q1 = (x1[kkk-1][i]-x1[kkk-1][j])**2 + (x2[kkk-1][i]-x2[kkk-1][j])**2 + a**2
            tao12_1 = ((x1[kkk-1][i]-x1[kkk-1][j])/q1)*(1-2*(x2[kkk-1][i]-x2[kkk-1][j])**2/q1)
            tao11_1 = -((x2[kkk-1][i]-x2[kkk-1][j])/q1)*(1+2*(x1[kkk-1][i]-x1[kkk-1][j])**2/q1)
            tao22_1 = ((x2[kkk-1][i]-x2[kkk-1][j])/q1)*(1-2*(x2[kkk-1][i]-x2[kkk-1][j])**2/q1)

            tao12_3 = ((x2[kkk-1][i]-x2[kkk-1][j])/q1)*(1-2*(x2[kkk-1][i]-x2[kkk-1][j])**2/q1)
            tao11_3 = ((x1[kkk-1][i]-x1[kkk-1][j])/q1)*(1-2*(x2[kkk-1][i]-x2[kkk-1][j])**2/q1)
            tao22_3 = ((x1[kkk-1][i]-x1[kkk-1][j])/q1)*(1+2*(x2[kkk-1][i]-x2[kkk-1][j])**2/q1)


            tao11_2 = coeff*(tao11_1+tao11_3)
            tao11_4 = coeff*(-tao11_1+tao11_3)
            tao12_2 = coeff*(tao12_1+tao12_3)
            tao12_4 = coeff*(-tao12_1+tao12_3)
            tao22_2 = coeff*(tao22_1+tao22_3)
            tao22_4 = coeff*(-tao22_1+tao22_3)

            if l1_old != 0 and l2_old/l1_old >= 0 and l2_old/l1_old < 1: 
                aaa = 0
                I11_12 += h11 * rho1_1[kkk-1][j] * length_jjj * tao12_1
                I12_12 += h12 * rho1_2[kkk-1][j] * length_jjj * tao12_2
                I11_11 += h11 * rho1_1[kkk-1][j] * length_jjj * tao11_1
                I12_11 += h12 * rho1_2[kkk-1][j] * length_jjj * tao11_2
                I11_22 += h11 * rho1_1[kkk-1][j] * length_jjj * tao22_1
                I12_22 += h12 * rho1_2[kkk-1][j] * length_jjj * tao22_2
                I21_12 += h21 * rho2_1[kkk-1][j] * length_jjj * tao12_1
                I22_12 += h22 * rho2_2[kkk-1][j] * length_jjj * tao12_2
                I21_11 += h21 * rho2_1[kkk-1][j] * length_jjj * tao11_1
                I22_11 += h22 * rho2_2[kkk-1][j] * length_jjj * tao11_2
                I21_22 += h21 * rho2_1[kkk-1][j] * length_jjj * tao22_1
                I22_22 += h22 * rho2_2[kkk-1][j] * length_jjj * tao22_2
                phi_ref1 = 0
                phi_ref2 = math.pi/ref
                beta11 = 1/(math.sin(phi_ref2-phi_ref1)) * beta110
                beta12 = 1/(math.sin(phi_ref2-phi_ref1)) * beta120
                beta13 = 1/(math.sin(phi_ref2-phi_ref1)) * beta130
                beta14 = 1/(math.sin(phi_ref2-phi_ref1)) * beta140 
                beta21 = 1/(math.sin(phi_ref2-phi_ref1)) * beta210
                beta22 = 1/(math.sin(phi_ref2-phi_ref1)) * beta220
                beta23 = 1/(math.sin(phi_ref2-phi_ref1)) * beta230
                beta24 = 1/(math.sin(phi_ref2-phi_ref1)) * beta240 
                LLL1 = np.array([1/(math.sin(phi_ref2-phi_ref1)) * beta110, 1/(math.sin(phi_ref2-phi_ref1)) * beta120])
                LLL2 = np.array([1/(math.sin(phi_ref2-phi_ref1)) * beta210, 1/(math.sin(phi_ref2-phi_ref1)) * beta220])
                if i != j:
                    f_R_1_1 += (2*rho1_1[kkk-1][j]*beta11*h11 + rho2_1[kkk-1][j]*beta21*h21 + rho1_2[kkk-1][j]*beta12*h12*coeff + rho1_4[kkk-1][j]*beta14*h14*coeff1 + rho2_2[kkk-1][j]*beta22*h22*coeff + rho2_4[kkk-1][j]*beta24*h24*coeff1)*np.log(RR/Ra)*length_jjj
                    f_R_2_1 += (rho1_1[kkk-1][j]*beta11*h11 + 2*rho2_1[kkk-1][j]*beta21*h21 + rho1_2[kkk-1][j]*beta12*h12*coeff + rho1_4[kkk-1][j]*beta14*h14*coeff1 + rho2_2[kkk-1][j]*beta22*h22*coeff + rho2_4[kkk-1][j]*beta24*h24*coeff1)*np.log(RR/Ra)*length_jjj
                    f_R_1_2 += (2*rho1_2[kkk-1][j]*beta12*h12 + rho2_2[kkk-1][j]*beta22*h22 + rho1_1[kkk-1][j]*beta11*h11*coeff + rho1_3[kkk-1][j]*beta13*h13*coeff + rho2_1[kkk-1][j]*beta21*h21*coeff + rho2_3[kkk-1][j]*beta23*h23*coeff)*np.log(RR/Ra)*length_jjj
                    f_R_2_2 += (rho1_2[kkk-1][j]*beta12*h12 + 2*rho2_2[kkk-1][j]*beta22*h22 + rho1_1[kkk-1][j]*beta11*h11*coeff + rho1_3[kkk-1][j]*beta13*h13*coeff + rho2_1[kkk-1][j]*beta21*h21*coeff + rho2_3[kkk-1][j]*beta23*h23*coeff)*np.log(RR/Ra)*length_jjj
    
                    f_th_1_1 += (2*rho1_1[kkk-1][j]*beta11*h11*costheta1*costheta1 + rho2_1[kkk-1][j]*beta21*h21*costheta1*costheta1 + rho1_2[kkk-1][j]*beta12*h12*costheta1*costheta2 + rho2_2[kkk-1][j]*beta22*h22*costheta1*costheta2 + rho1_3[kkk-1][j]*beta13*h13*costheta1*costheta3 + rho2_3[kkk-1][j]*beta23*h23*costheta1*costheta3 + rho1_4[kkk-1][j]*beta14*h14*costheta1*costheta4 + rho2_4[kkk-1][j]*beta24*h24*costheta1*costheta4)*length_jjj
                    f_th_2_1 += (rho1_1[kkk-1][j]*beta11*h11*costheta1*costheta1 + 2*rho2_1[kkk-1][j]*beta21*h21*costheta1*costheta1 + rho1_2[kkk-1][j]*beta12*h12*costheta1*costheta2 + rho2_2[kkk-1][j]*beta22*h22*costheta1*costheta2 + rho1_3[kkk-1][j]*beta13*h13*costheta1*costheta3 + rho2_3[kkk-1][j]*beta23*h23*costheta1*costheta3 + rho1_4[kkk-1][j]*beta14*h14*costheta1*costheta4 + rho2_4[kkk-1][j]*beta24*h24*costheta1*costheta4)*length_jjj
                    f_th_1_2 += (rho1_1[kkk-1][j]*beta11*h11*costheta1*costheta2 + rho2_1[kkk-1][j]*beta21*h21*costheta1*costheta2 + 2*rho1_2[kkk-1][j]*beta12*h12*costheta2*costheta2 + rho2_2[kkk-1][j]*beta22*h22*costheta2*costheta2 + rho1_3[kkk-1][j]*beta13*h13*costheta2*costheta3 + rho2_3[kkk-1][j]*beta23*h23*costheta2*costheta3 + rho1_4[kkk-1][j]*beta14*h14*costheta2*costheta4 + rho2_4[kkk-1][j]*beta24*h24*costheta2*costheta4)*length_jjj
                    f_th_2_2 += (rho1_1[kkk-1][j]*beta11*h11*costheta1*costheta2 + rho2_1[kkk-1][j]*beta21*h21*costheta1*costheta2 + rho1_2[kkk-1][j]*beta12*h12*costheta2*costheta2 + 2*rho2_2[kkk-1][j]*beta22*h22*costheta2*costheta2 + rho1_3[kkk-1][j]*beta13*h13*costheta2*costheta3 + rho2_3[kkk-1][j]*beta23*h23*costheta2*costheta3 + rho1_4[kkk-1][j]*beta14*h14*costheta2*costheta4 + rho2_4[kkk-1][j]*beta24*h24*costheta2*costheta4)*length_jjj


            if l1_old != 0 and l2_old/l1_old >= 1:
                aaa = 1
                I12_12 += h12 * rho1_2[kkk-1][j] * length_jjj * tao12_2
                I13_12 += h13 * rho1_3[kkk-1][j] * length_jjj * tao12_3
                I12_11 += h12 * rho1_2[kkk-1][j] * length_jjj * tao11_2
                I13_11 += h13 * rho1_3[kkk-1][j] * length_jjj * tao11_3 
                I12_22 += h12 * rho1_2[kkk-1][j] * length_jjj * tao22_2
                I13_22 += h13 * rho1_3[kkk-1][j] * length_jjj * tao22_3
                I22_12 += h22 * rho2_2[kkk-1][j] * length_jjj * tao12_2
                I23_12 += h23 * rho2_3[kkk-1][j] * length_jjj * tao12_3
                I22_11 += h22 * rho2_2[kkk-1][j] * length_jjj * tao11_2
                I23_11 += h23 * rho2_3[kkk-1][j] * length_jjj * tao11_3 
                I22_22 += h22 * rho2_2[kkk-1][j] * length_jjj * tao22_2
                I23_22 += h23 * rho2_3[kkk-1][j] * length_jjj * tao22_3

                phi_ref1 = math.pi/ref
                phi_ref2 = math.pi/ref * 2
                beta11 = 1/(math.sin(phi_ref2-phi_ref1)) * beta110
                beta12 = 1/(math.sin(phi_ref2-phi_ref1)) * beta120
                beta13 = 1/(math.sin(phi_ref2-phi_ref1)) * beta130
                beta14 = 1/(math.sin(phi_ref2-phi_ref1)) * beta140 
                beta21 = 1/(math.sin(phi_ref2-phi_ref1)) * beta210
                beta22 = 1/(math.sin(phi_ref2-phi_ref1)) * beta220
                beta23 = 1/(math.sin(phi_ref2-phi_ref1)) * beta230
                beta24 = 1/(math.sin(phi_ref2-phi_ref1)) * beta240 
                LLL1 = np.array([1/(math.sin(phi_ref2-phi_ref1)) * beta120, 1/(math.sin(phi_ref2-phi_ref1)) * beta130])
                LLL2 = np.array([1/(math.sin(phi_ref2-phi_ref1)) * beta220, 1/(math.sin(phi_ref2-phi_ref1)) * beta230])
                if i != j:
                    f_R_1_2 += (2*rho1_2[kkk-1][j]*beta12*h12 + rho2_2[kkk-1][j]*beta22*h22 + rho1_1[kkk-1][j]*beta11*h11*coeff + rho1_3[kkk-1][j]*beta13*h13*coeff + rho2_1[kkk-1][j]*beta21*h21*coeff + rho2_3[kkk-1][j]*beta23*h23*coeff)*np.log(RR/Ra)*length_jjj
                    f_R_2_2 += (rho1_2[kkk-1][j]*beta12*h12 + 2*rho2_2[kkk-1][j]*beta22*h22 + rho1_1[kkk-1][j]*beta11*h11*coeff + rho1_3[kkk-1][j]*beta13*h13*coeff + rho2_1[kkk-1][j]*beta21*h21*coeff + rho2_3[kkk-1][j]*beta23*h23*coeff)*np.log(RR/Ra)*length_jjj
                    f_R_1_3 += (2*rho1_3[kkk-1][j]*beta13*h13 + rho2_3[kkk-1][j]*beta23*h23 + rho1_2[kkk-1][j]*beta12*h12*coeff + rho1_4[kkk-1][j]*beta14*h14*coeff + rho2_2[kkk-1][j]*beta22*h22*coeff + rho2_4[kkk-1][j]*beta24*h24*coeff)*np.log(RR/Ra)*length_jjj
                    f_R_2_3 += (rho1_3[kkk-1][j]*beta13*h13 + 2*rho2_3[kkk-1][j]*beta23*h23 + rho1_2[kkk-1][j]*beta12*h12*coeff + rho1_4[kkk-1][j]*beta14*h14*coeff + rho2_2[kkk-1][j]*beta22*h22*coeff + rho2_4[kkk-1][j]*beta24*h24*coeff)*np.log(RR/Ra)*length_jjj

                    f_th_1_2 += (rho1_1[kkk-1][j]*beta11*h11*costheta1*costheta2 + rho2_1[kkk-1][j]*beta21*h21*costheta1*costheta2 + 2*rho1_2[kkk-1][j]*beta12*h12*costheta2*costheta2 + rho2_2[kkk-1][j]*beta22*h22*costheta2*costheta2 + rho1_3[kkk-1][j]*beta13*h13*costheta2*costheta3 + rho2_3[kkk-1][j]*beta23*h23*costheta2*costheta3 + rho1_4[kkk-1][j]*beta14*h14*costheta2*costheta4 + rho2_4[kkk-1][j]*beta24*h24*costheta2*costheta4)*length_jjj
                    f_th_2_2 += (rho1_1[kkk-1][j]*beta11*h11*costheta1*costheta2 + rho2_1[kkk-1][j]*beta21*h21*costheta1*costheta2 + rho1_2[kkk-1][j]*beta12*h12*costheta2*costheta2 + 2*rho2_2[kkk-1][j]*beta22*h22*costheta2*costheta2 + rho1_3[kkk-1][j]*beta13*h13*costheta2*costheta3 + rho2_3[kkk-1][j]*beta23*h23*costheta2*costheta3 + rho1_4[kkk-1][j]*beta14*h14*costheta2*costheta4 + rho2_4[kkk-1][j]*beta24*h24*costheta2*costheta4)*length_jjj
                    f_th_1_3 += (rho1_1[kkk-1][j]*beta11*h11*costheta3*costheta1 + rho2_1[kkk-1][j]*beta21*h21*costheta3*costheta1 + rho1_2[kkk-1][j]*beta12*h12*costheta3*costheta2 + rho2_2[kkk-1][j]*beta22*h22*costheta3*costheta2 + 2*rho1_3[kkk-1][j]*beta13*h13*costheta3*costheta3 + rho2_3[kkk-1][j]*beta23*h23*costheta3*costheta3 + rho1_4[kkk-1][j]*beta14*h14*costheta3*costheta4 + rho2_4[kkk-1][j]*beta24*h24*costheta3*costheta4)*length_jjj
                    f_th_2_3 += (rho1_1[kkk-1][j]*beta11*h11*costheta3*costheta1 + rho2_1[kkk-1][j]*beta21*h21*costheta3*costheta1 + rho1_2[kkk-1][j]*beta12*h12*costheta3*costheta2 + rho2_2[kkk-1][j]*beta22*h22*costheta3*costheta2 + rho1_3[kkk-1][j]*beta13*h13*costheta3*costheta3 + 2*rho2_3[kkk-1][j]*beta23*h23*costheta3*costheta3 + rho1_4[kkk-1][j]*beta14*h14*costheta3*costheta4 + rho2_4[kkk-1][j]*beta24*h24*costheta3*costheta4)*length_jjj



            if l1_old == 0 or (l1_old != 0 and l2_old/l1_old < -1):              
                aaa = 2
                I13_12 += h13 * rho1_3[kkk-1][j] * length_jjj * tao12_3
                I14_12 += h14 * rho1_4[kkk-1][j] * length_jjj * tao12_4
                I13_11 += h13 * rho1_3[kkk-1][j] * length_jjj * tao11_3
                I14_11 += h14 * rho1_4[kkk-1][j] * length_jjj * tao11_4    
                I13_22 += h13 * rho1_3[kkk-1][j] * length_jjj * tao22_3
                I14_22 += h14 * rho1_4[kkk-1][j] * length_jjj * tao22_4              
                I23_12 += h23 * rho2_3[kkk-1][j] * length_jjj * tao12_3
                I24_12 += h24 * rho2_4[kkk-1][j] * length_jjj * tao12_4
                I23_11 += h23 * rho2_3[kkk-1][j] * length_jjj * tao11_3
                I24_11 += h24 * rho2_4[kkk-1][j] * length_jjj * tao11_4    
                I23_22 += h23 * rho2_3[kkk-1][j] * length_jjj * tao22_3
                I24_22 += h24 * rho2_4[kkk-1][j] * length_jjj * tao22_4              

                phi_ref1 = math.pi/ref * 2
                phi_ref2 = math.pi/ref * 3   
                beta11 = 1/(math.sin(phi_ref2-phi_ref1)) * beta110
                beta12 = 1/(math.sin(phi_ref2-phi_ref1)) * beta120
                beta13 = 1/(math.sin(phi_ref2-phi_ref1)) * beta130
                beta14 = 1/(math.sin(phi_ref2-phi_ref1)) * beta140 
                beta21 = 1/(math.sin(phi_ref2-phi_ref1)) * beta210
                beta22 = 1/(math.sin(phi_ref2-phi_ref1)) * beta220
                beta23 = 1/(math.sin(phi_ref2-phi_ref1)) * beta230
                beta24 = 1/(math.sin(phi_ref2-phi_ref1)) * beta240 
                LLL1 = np.array([1/(math.sin(phi_ref2-phi_ref1)) * beta130, 1/(math.sin(phi_ref2-phi_ref1)) * beta140])
                LLL2 = np.array([1/(math.sin(phi_ref2-phi_ref1)) * beta230, 1/(math.sin(phi_ref2-phi_ref1)) * beta240])
                if i != j:
                    f_R_1_3 += (2*rho1_3[kkk-1][j]*beta13*h13 + rho2_3[kkk-1][j]*beta23*h23 + rho1_2[kkk-1][j]*beta12*h12*coeff + rho1_4[kkk-1][j]*beta14*h14*coeff + rho2_2[kkk-1][j]*beta22*h22*coeff + rho2_4[kkk-1][j]*beta24*h24*coeff)*np.log(RR/Ra)*length_jjj
                    f_R_2_3 += (rho1_3[kkk-1][j]*beta13*h13 + 2*rho2_3[kkk-1][j]*beta23*h23 + rho1_2[kkk-1][j]*beta12*h12*coeff + rho1_4[kkk-1][j]*beta14*h14*coeff + rho2_2[kkk-1][j]*beta22*h22*coeff + rho2_4[kkk-1][j]*beta24*h24*coeff)*np.log(RR/Ra)*length_jjj
                    f_R_1_4 += (2*rho1_4[kkk-1][j]*beta14*h14 + rho2_4[kkk-1][j]*beta24*h24 + rho1_1[kkk-1][j]*beta11*h11*coeff1 + rho1_3[kkk-1][j]*beta13*h13*coeff + rho2_1[kkk-1][j]*beta21*h21*coeff1 + rho2_3[kkk-1][j]*beta23*h23*coeff)*np.log(RR/Ra)*length_jjj
                    f_R_2_4 += (rho1_4[kkk-1][j]*beta14*h14 + 2*rho2_4[kkk-1][j]*beta24*h24 + rho1_1[kkk-1][j]*beta11*h11*coeff1 + rho1_3[kkk-1][j]*beta13*h13*coeff + rho2_1[kkk-1][j]*beta21*h21*coeff1 + rho2_3[kkk-1][j]*beta23*h23*coeff)*np.log(RR/Ra)*length_jjj

                    f_th_1_3 += (rho1_1[kkk-1][j]*beta11*h11*costheta3*costheta1 + rho2_1[kkk-1][j]*beta21*h21*costheta3*costheta1 + rho1_2[kkk-1][j]*beta12*h12*costheta3*costheta2 + rho2_2[kkk-1][j]*beta22*h22*costheta3*costheta2 + 2*rho1_3[kkk-1][j]*beta13*h13*costheta3*costheta3 + rho2_3[kkk-1][j]*beta23*h23*costheta3*costheta3 + rho1_4[kkk-1][j]*beta14*h14*costheta3*costheta4 + rho2_4[kkk-1][j]*beta24*h24*costheta3*costheta4)*length_jjj
                    f_th_2_3 += (rho1_1[kkk-1][j]*beta11*h11*costheta3*costheta1 + rho2_1[kkk-1][j]*beta21*h21*costheta3*costheta1 + rho1_2[kkk-1][j]*beta12*h12*costheta3*costheta2 + rho2_2[kkk-1][j]*beta22*h22*costheta3*costheta2 + rho1_3[kkk-1][j]*beta13*h13*costheta3*costheta3 + 2*rho2_3[kkk-1][j]*beta23*h23*costheta3*costheta3 + rho1_4[kkk-1][j]*beta14*h14*costheta3*costheta4 + rho2_4[kkk-1][j]*beta24*h24*costheta3*costheta4)*length_jjj
                    f_th_1_4 += (rho1_1[kkk-1][j]*beta11*h11*costheta4*costheta1 + rho2_1[kkk-1][j]*beta21*h21*costheta4*costheta1 + rho1_2[kkk-1][j]*beta12*h12*costheta4*costheta2 + rho2_2[kkk-1][j]*beta22*h22*costheta4*costheta2 + rho1_3[kkk-1][j]*beta13*h13*costheta4*costheta3 + rho2_3[kkk-1][j]*beta23*h23*costheta4*costheta3 + 2*rho1_4[kkk-1][j]*beta14*h14*costheta4*costheta4 + rho2_4[kkk-1][j]*beta24*h24*costheta4*costheta4)*length_jjj
                    f_th_2_4 += (rho1_1[kkk-1][j]*beta11*h11*costheta4*costheta1 + rho2_1[kkk-1][j]*beta21*h21*costheta4*costheta1 + rho1_2[kkk-1][j]*beta12*h12*costheta4*costheta2 + rho2_2[kkk-1][j]*beta22*h22*costheta4*costheta2 + rho1_3[kkk-1][j]*beta13*h13*costheta4*costheta3 + rho2_3[kkk-1][j]*beta23*h23*costheta4*costheta3 + rho1_4[kkk-1][j]*beta14*h14*costheta4*costheta4 + 2*rho2_4[kkk-1][j]*beta24*h24*costheta4*costheta4)*length_jjj


            if l1_old != 0 and l2_old/l1_old >= -1 and l2_old/l1_old < 0:                
                aaa = 3
                I14_12 +=  -h14 * rho1_4[kkk-1][j] * length_jjj * tao12_4
                I11_12 +=  -h11 * rho1_1[kkk-1][j] * length_jjj * tao12_1   
                I14_11 +=  -h14 * rho1_4[kkk-1][j] * length_jjj * tao11_4
                I11_11 +=  -h11 * rho1_1[kkk-1][j] * length_jjj * tao11_1                   
                I14_22 +=  -h14 * rho1_4[kkk-1][j] * length_jjj * tao22_4
                I11_22 +=  -h11 * rho1_1[kkk-1][j] * length_jjj * tao22_1         
                I24_12 +=  -h24 * rho2_4[kkk-1][j] * length_jjj * tao12_4
                I21_12 +=  -h21 * rho2_1[kkk-1][j] * length_jjj * tao12_1   
                I24_11 +=  -h24 * rho2_4[kkk-1][j] * length_jjj * tao11_4
                I21_11 +=  -h21 * rho2_1[kkk-1][j] * length_jjj * tao11_1                   
                I24_22 +=  -h24 * rho2_4[kkk-1][j] * length_jjj * tao22_4
                I21_22 +=  -h21 * rho2_1[kkk-1][j] * length_jjj * tao22_1         

                phi_ref1 = math.pi/ref * 3
                phi_ref2 = math.pi/ref * 0
                beta11 = 1/(math.sin(phi_ref2-phi_ref1)) * beta110
                beta12 = 1/(math.sin(phi_ref2-phi_ref1)) * beta120
                beta13 = 1/(math.sin(phi_ref2-phi_ref1)) * beta130
                beta14 = 1/(math.sin(phi_ref2-phi_ref1)) * beta140 
                beta21 = 1/(math.sin(phi_ref2-phi_ref1)) * beta210
                beta22 = 1/(math.sin(phi_ref2-phi_ref1)) * beta220
                beta23 = 1/(math.sin(phi_ref2-phi_ref1)) * beta230
                beta24 = 1/(math.sin(phi_ref2-phi_ref1)) * beta240 
                LLL1 = np.array([1/(math.sin(phi_ref2-phi_ref1)) * beta140, 1/(math.sin(phi_ref2-phi_ref1)) * beta110])
                LLL2 = np.array([1/(math.sin(phi_ref2-phi_ref1)) * beta240, 1/(math.sin(phi_ref2-phi_ref1)) * beta210])                                
                if i != j:
                    f_R_1_4 += (2*rho1_4[kkk-1][j]*beta14*h14 + rho2_4[kkk-1][j]*beta24*h24 + rho1_1[kkk-1][j]*beta11*h11*coeff1 + rho1_3[kkk-1][j]*beta13*h13*coeff + rho2_1[kkk-1][j]*beta21*h21*coeff1 + rho2_3[kkk-1][j]*beta23*h23*coeff)*np.log(RR/Ra)*length_jjj
                    f_R_2_4 += (rho1_4[kkk-1][j]*beta14*h14 + 2*rho2_4[kkk-1][j]*beta24*h24 + rho1_1[kkk-1][j]*beta11*h11*coeff1 + rho1_3[kkk-1][j]*beta13*h13*coeff + rho2_1[kkk-1][j]*beta21*h21*coeff1 + rho2_3[kkk-1][j]*beta23*h23*coeff)*np.log(RR/Ra)*length_jjj
                    f_R_1_1 += (2*rho1_1[kkk-1][j]*beta11*h11 + rho2_1[kkk-1][j]*beta21*h21 + rho1_2[kkk-1][j]*beta12*h12*coeff + rho1_4[kkk-1][j]*beta14*h14*coeff1 + rho2_2[kkk-1][j]*beta22*h22*coeff + rho2_4[kkk-1][j]*beta24*h24*coeff1)*np.log(RR/Ra)*length_jjj
                    f_R_2_1 += (rho1_1[kkk-1][j]*beta11*h11 + 2*rho2_1[kkk-1][j]*beta21*h21 + rho1_2[kkk-1][j]*beta12*h12*coeff + rho1_4[kkk-1][j]*beta14*h14*coeff1 + rho2_2[kkk-1][j]*beta22*h22*coeff + rho2_4[kkk-1][j]*beta24*h24*coeff1)*np.log(RR/Ra)*length_jjj

                    f_th_1_1 += (2*rho1_1[kkk-1][j]*beta11*h11*costheta1*costheta1 + rho2_1[kkk-1][j]*beta21*h21*costheta1*costheta1 + rho1_2[kkk-1][j]*beta12*h12*costheta1*costheta2 + rho2_2[kkk-1][j]*beta22*h22*costheta1*costheta2 + rho1_3[kkk-1][j]*beta13*h13*costheta1*costheta3 + rho2_3[kkk-1][j]*beta23*h23*costheta1*costheta3 + rho1_4[kkk-1][j]*beta14*h14*costheta1*costheta4 + rho2_4[kkk-1][j]*beta24*h24*costheta1*costheta4)*length_jjj
                    f_th_2_1 += (rho1_1[kkk-1][j]*beta11*h11*costheta1*costheta1 + 2*rho2_1[kkk-1][j]*beta21*h21*costheta1*costheta1 + rho1_2[kkk-1][j]*beta12*h12*costheta1*costheta2 + rho2_2[kkk-1][j]*beta22*h22*costheta1*costheta2 + rho1_3[kkk-1][j]*beta13*h13*costheta1*costheta3 + rho2_3[kkk-1][j]*beta23*h23*costheta1*costheta3 + rho1_4[kkk-1][j]*beta14*h14*costheta1*costheta4 + rho2_4[kkk-1][j]*beta24*h24*costheta1*costheta4)*length_jjj
                    f_th_1_4 += (rho1_1[kkk-1][j]*beta11*h11*costheta4*costheta1 + rho2_1[kkk-1][j]*beta21*h21*costheta4*costheta1 + rho1_2[kkk-1][j]*beta12*h12*costheta4*costheta2 + rho2_2[kkk-1][j]*beta22*h22*costheta4*costheta2 + rho1_3[kkk-1][j]*beta13*h13*costheta4*costheta3 + rho2_3[kkk-1][j]*beta23*h23*costheta4*costheta3 + 2*rho1_4[kkk-1][j]*beta14*h14*costheta4*costheta4 + rho2_4[kkk-1][j]*beta24*h24*costheta4*costheta4)*length_jjj
                    f_th_2_4 += (rho1_1[kkk-1][j]*beta11*h11*costheta4*costheta1 + rho2_1[kkk-1][j]*beta21*h21*costheta4*costheta1 + rho1_2[kkk-1][j]*beta12*h12*costheta4*costheta2 + rho2_2[kkk-1][j]*beta22*h22*costheta4*costheta2 + rho1_3[kkk-1][j]*beta13*h13*costheta4*costheta3 + rho2_3[kkk-1][j]*beta23*h23*costheta4*costheta3 + rho1_4[kkk-1][j]*beta14*h14*costheta4*costheta4 + 2*rho2_4[kkk-1][j]*beta24*h24*costheta4*costheta4)*length_jjj

            
        sigma11_1 = G*(beta11 * I11_11 + beta12 * I12_11 + beta13 * I13_11 + beta14 * I14_11)
        sigma12_1 = G*(beta11 * I11_12 + beta12 * I12_12 + beta13 * I13_12 + beta14 * I14_12) #+ tao_ext
        sigma22_1 = G*(beta11 * I11_22 + beta12 * I12_22 + beta13 * I13_22 + beta14 * I14_22)   
        sigma11_2 = G*(beta21 * I21_11 + beta22 * I22_11 + beta23 * I23_11 + beta24 * I24_11)
        sigma12_2 = G*(beta21 * I21_12 + beta22 * I22_12 + beta23 * I23_12 + beta24 * I24_12) #+ tao_ext
        sigma22_2 = G*(beta21 * I21_22 + beta22 * I22_22 + beta23 * I23_22 + beta24 * I24_22)   

        sigma11 = sigma11_1 + sigma11_2
        sigma12 = sigma12_1 + sigma12_2
        sigma22 = sigma22_1 + sigma22_2

        tao_k = (sigma22-sigma11)/2 * math.sin(2*phi_ref1) + sigma12 * math.cos(2*phi_ref1)
        tao_k1 = (sigma22-sigma11)/2 * math.sin(2*phi_ref2) + sigma12 * math.cos(2*phi_ref2)
        

        tao1_k = (sigma22_1-sigma11_1)/2 * math.sin(2*phi_ref1) + sigma12_1 * math.cos(2*phi_ref1)
        tao1_k1 = (sigma22_1-sigma11_1)/2 * math.sin(2*phi_ref2) + sigma12_1 * math.cos(2*phi_ref2)
        tao1 = np.array([-tao1_k, -tao1_k1])        
        tao2_k = (sigma22_2-sigma11_2)/2 * math.sin(2*phi_ref1) + sigma12_2 * math.cos(2*phi_ref1)
        tao2_k1 = (sigma22_2-sigma11_2)/2 * math.sin(2*phi_ref2) + sigma12_2 * math.cos(2*phi_ref2)
        tao2 = np.array([-tao2_k, -tao2_k1])        

        kap = (dev1_1*dev2_2-dev2_1*dev1_2)/(dev1_1**2+dev1_2**2)**1.5


        if l1_old != 0 and l2_old/l1_old >= 0 and l2_old/l1_old < 1:
            gama_rho11 = gamma1*np.sign(h11*rho1_1[kkk-1][i] + h21*rho2_1[kkk-1][i])*sca*h11
            gama_rho12 = gamma2*np.sign(h12*rho1_2[kkk-1][i] + h22*rho2_2[kkk-1][i])*sca*h12
            gama_rho21 = gamma1*np.sign(h11*rho1_1[kkk-1][i] + h21*rho2_1[kkk-1][i])*sca*h21
            gama_rho22 = gamma2*np.sign(h12*rho1_2[kkk-1][i] + h22*rho2_2[kkk-1][i])*sca*h22
            gama_rho13 = 0
            gama_rho14 = 0
            gama_rho23 = 0
            gama_rho24 = 0
            f_rho_1_1 = stablize(G_r*(beta11*h11*(f_R_1_1 + f_th_1_1)) + LL_sca*lamda_new_value * (2*h11*(h11*rho1_1[kkk-1][i] + h21*rho2_1[kkk-1][i])-2*h11*(h12*rho1_2[kkk-1][i]+h22*rho2_2[kkk-1][i])*coeff)/length_kkk - gama_rho11)
            f_rho_2_1 = stablize(G_r*(beta21*h21*(f_R_2_1 + f_th_2_1)) + LL_sca*lamda_new_value * (2*h21*(h11*rho1_1[kkk-1][i] + h21*rho2_1[kkk-1][i])-2*h21*(h12*rho1_2[kkk-1][i]+h22*rho2_2[kkk-1][i])*coeff)/length_kkk - gama_rho21)
            f_rho_1_2 = stablize(G_r*(beta12*h12*(f_R_1_2 + f_th_1_2)) + LL_sca*lamda_new_value * (2*h12*(h12*rho1_2[kkk-1][i] + h22*rho2_2[kkk-1][i])-2*h12*(h11*rho1_1[kkk-1][i]+h21*rho2_1[kkk-1][i])*coeff)/length_kkk - gama_rho12)
            f_rho_2_2 = stablize(G_r*(beta22*h22*(f_R_2_2 + f_th_2_2)) + LL_sca*lamda_new_value * (2*h22*(h12*rho1_2[kkk-1][i] + h22*rho2_2[kkk-1][i])-2*h22*(h11*rho1_1[kkk-1][i]+h21*rho2_1[kkk-1][i])*coeff)/length_kkk - gama_rho22)
            f_rho_1_3 = 0
            f_rho_2_3 = 0
            f_rho_1_4 = 0
            f_rho_2_4 = 0
            rho11_new_value = rho1_1[kkk-1][i] + delt*M_rho11*f_rho_1_1
            rho21_new_value = rho2_1[kkk-1][i] + delt*M_rho21*f_rho_2_1
            rho12_new_value = rho1_2[kkk-1][i] + delt*M_rho12*f_rho_1_2
            rho22_new_value = rho2_2[kkk-1][i] + delt*M_rho22*f_rho_2_2
            rho13_new_value = rho1_3[kkk-1][i] + delt*M_rho11*f_rho_1_3
            rho23_new_value = rho2_3[kkk-1][i] + delt*M_rho21*f_rho_2_3
            rho14_new_value = rho1_4[kkk-1][i] + delt*M_rho12*f_rho_1_4
            rho24_new_value = rho2_4[kkk-1][i] + delt*M_rho22*f_rho_2_4
                        
            l1_new = - (h12 * rho12_new_value + h22 * rho22_new_value) * np.cos(0) + (h11 * rho11_new_value + h21 * rho21_new_value) * np.cos(np.pi/4)
            l2_new = - (h12 * rho12_new_value + h22 * rho22_new_value) * np.sin(0) + (h11 * rho11_new_value + h21 * rho21_new_value) * np.sin(np.pi/4)
            f_g1 = -(M_11*tao1_k*beta11*h11*rho11_new_value + M_21*tao2_k*beta21*h21*rho21_new_value) - (stablize_l(np.sign(h11*rho11_new_value + h21*rho21_new_value)/(h12*rho12_new_value + h22*rho22_new_value))*(beta12*h12*rho12_new_value + beta22*h22*rho22_new_value)*(M_11*tao1_k1*abs(h11*rho11_new_value)+M_21*tao2_k1*abs(h21*rho21_new_value))  )
            f_g2 = -(M_12*tao1_k1*beta12*h12*rho12_new_value + M_22*tao2_k1*beta22*h22*rho22_new_value) - (stablize_l(np.sign(h12*rho12_new_value + h22*rho22_new_value)/(h11*rho11_new_value + h21*rho21_new_value))*(beta11*h11*rho11_new_value + beta21*h21*rho21_new_value)*(M_12*tao1_k*abs(h12*rho12_new_value)+M_22*tao2_k*abs(h22*rho22_new_value))  )
            f_g1_r = f_g1 * np.cos(np.pi/4*0) + f_g2 * np.cos(np.pi/4*1)
            f_g2_r = f_g1 * np.sin(np.pi/4*0) + f_g2 * np.sin(np.pi/4)
            l_1 = -(h12 * rho12_new_value + h22 * rho22_new_value)
            l_2 = (h11 * rho11_new_value + h21 * rho21_new_value)

        if l1_old != 0 and l2_old/l1_old >= 1:
            gama_rho11 = 0
            gama_rho21 = 0
            gama_rho12 = gamma2*np.sign(h12*rho1_2[kkk-1][i] + h22*rho2_2[kkk-1][i])*sca*h12
            gama_rho22 = gamma2*np.sign(h12*rho1_2[kkk-1][i] + h22*rho2_2[kkk-1][i])*sca*h22
            gama_rho13 = gamma3*np.sign(h13*rho1_3[kkk-1][i] + h23*rho2_3[kkk-1][i])*sca*h13
            gama_rho23 = gamma3*np.sign(h13*rho1_3[kkk-1][i] + h23*rho2_3[kkk-1][i])*sca*h23
            gama_rho14 = 0
            gama_rho24 = 0
            f_rho_1_1 = 0
            f_rho_2_1 = 0
            f_rho_1_2 = stablize(G_r*(beta12*h12*(f_R_1_2 + f_th_1_2)) + LL_sca*lamda_new_value * (2*h12*(h12*rho1_2[kkk-1][i] + h22*rho2_2[kkk-1][i])-2*h12*(h13*rho1_3[kkk-1][i]+h23*rho2_3[kkk-1][i])*coeff)/length_kkk - gama_rho12)
            f_rho_2_2 = stablize(G_r*(beta22*h22*(f_R_2_2 + f_th_2_2)) + LL_sca*lamda_new_value * (2*h22*(h12*rho1_2[kkk-1][i] + h22*rho2_2[kkk-1][i])-2*h22*(h13*rho1_3[kkk-1][i]+h23*rho2_3[kkk-1][i])*coeff)/length_kkk - gama_rho22)
            f_rho_1_3 = stablize(G_r*(beta13*h13*(f_R_1_3 + f_th_1_3)) + LL_sca*lamda_new_value * (2*h13*(h13*rho1_3[kkk-1][i] + h23*rho2_3[kkk-1][i])-2*h13*(h12*rho1_2[kkk-1][i]+h22*rho2_2[kkk-1][i])*coeff)/length_kkk - gama_rho13)
            f_rho_2_3 = stablize(G_r*(beta23*h23*(f_R_2_3 + f_th_2_3)) + LL_sca*lamda_new_value * (2*h23*(h13*rho1_3[kkk-1][i] + h23*rho2_3[kkk-1][i])-2*h23*(h12*rho1_2[kkk-1][i]+h22*rho2_2[kkk-1][i])*coeff)/length_kkk - gama_rho23)
            f_rho_1_4 = 0
            f_rho_2_4 = 0
            rho11_new_value = rho1_1[kkk-1][i] + delt*M_rho11*f_rho_1_1
            rho21_new_value = rho2_1[kkk-1][i] + delt*M_rho21*f_rho_2_1
            rho12_new_value = rho1_2[kkk-1][i] + delt*M_rho12*f_rho_1_2
            rho22_new_value = rho2_2[kkk-1][i] + delt*M_rho22*f_rho_2_2
            rho13_new_value = rho1_3[kkk-1][i] + delt*M_rho11*f_rho_1_3
            rho23_new_value = rho2_3[kkk-1][i] + delt*M_rho21*f_rho_2_3
            rho14_new_value = rho1_4[kkk-1][i] + delt*M_rho12*f_rho_1_4
            rho24_new_value = rho2_4[kkk-1][i] + delt*M_rho22*f_rho_2_4

            l1_new = - (h13 * rho13_new_value + h23 * rho23_new_value) * np.cos(np.pi/4) + (h12 * rho12_new_value + h22 * rho22_new_value) * np.cos(np.pi/2)
            l2_new = - (h13 * rho13_new_value + h23 * rho23_new_value) * np.sin(np.pi/4) + (h12 * rho12_new_value + h22 * rho22_new_value) * np.sin(np.pi/2)
            f_g1 = -(M_12*tao1_k*beta12*h12*rho12_new_value + M_22*tao2_k* beta22*h22*rho22_new_value) - (stablize_l(np.sign(h12*rho12_new_value + h22*rho22_new_value)/(h13*rho13_new_value + h23*rho23_new_value))*(beta13*h13*rho13_new_value + beta23*h23*rho23_new_value)*(M_12*tao1_k1*abs(h12*rho12_new_value)+M_22*tao2_k1*abs(h22*rho22_new_value))  )
            f_g2 = -(M_11*tao1_k1*beta13*h13*rho13_new_value + M_21*tao2_k1*beta23*h23*rho23_new_value) - (stablize_l(np.sign(h13*rho13_new_value + h23*rho23_new_value)/(h12*rho12_new_value + h22*rho22_new_value))*(beta12*h12*rho12_new_value + beta22*h22*rho22_new_value)*(M_11*tao1_k*abs(h13*rho13_new_value)+M_21*tao2_k*abs(h23*rho23_new_value))  )
            f_g1_r = f_g1 * np.cos(1*np.pi/4) + f_g2 * np.cos(2*np.pi/4)
            f_g2_r = f_g1 * np.sin(1*np.pi/4) + f_g2 * np.sin(2*np.pi/4)
            l_1 = -(h13 * rho13_new_value + h23 * rho23_new_value)
            l_2 = (h12 * rho12_new_value + h22 * rho22_new_value)

        if l1_old == 0 or (l1_old != 0 and l2_old/l1_old < -1): 
            gama_rho11 = 0
            gama_rho21 = 0
            gama_rho12 = 0
            gama_rho22 = 0
            gama_rho13 = gamma3*np.sign(h13*rho1_3[kkk-1][i] + h23*rho2_3[kkk-1][i])*sca*h13
            gama_rho23 = gamma3*np.sign(h13*rho1_3[kkk-1][i] + h23*rho2_3[kkk-1][i])*sca*h23
            gama_rho14 = gamma4*np.sign(h14*rho1_4[kkk-1][i] + h24*rho2_4[kkk-1][i])*sca*h14
            gama_rho24 = gamma4*np.sign(h14*rho1_4[kkk-1][i] + h24*rho2_4[kkk-1][i])*sca*h24
            f_rho_1_1 = 0
            f_rho_2_1 = 0
            f_rho_1_2 = 0
            f_rho_2_2 = 0
            f_rho_1_3 = stablize(G_r*(beta13*h13*(f_R_1_3 + f_th_1_3)) + LL_sca*lamda_new_value * (2*h13*(h13*rho1_3[kkk-1][i] + h23*rho2_3[kkk-1][i])-2*h13*(h14*rho1_4[kkk-1][i]+h24*rho2_4[kkk-1][i])*coeff)/length_kkk - gama_rho13)
            f_rho_2_3 = stablize(G_r*(beta23*h23*(f_R_2_3 + f_th_2_3)) + LL_sca*lamda_new_value * (2*h23*(h13*rho1_3[kkk-1][i] + h23*rho2_3[kkk-1][i])-2*h23*(h14*rho1_4[kkk-1][i]+h24*rho2_4[kkk-1][i])*coeff)/length_kkk - gama_rho23)
            f_rho_1_4 = stablize(G_r*(beta14*h14*(f_R_1_4 + f_th_1_4)) + LL_sca*lamda_new_value * (2*h14*(h14*rho1_4[kkk-1][i] + h24*rho2_4[kkk-1][i])-2*h14*(h13*rho1_3[kkk-1][i]+h23*rho2_3[kkk-1][i])*coeff)/length_kkk - gama_rho14)
            f_rho_2_4 = stablize(G_r*(beta24*h24*(f_R_2_4 + f_th_2_4)) + LL_sca*lamda_new_value * (2*h24*(h14*rho1_4[kkk-1][i] + h24*rho2_4[kkk-1][i])-2*h24*(h13*rho1_3[kkk-1][i]+h23*rho2_3[kkk-1][i])*coeff)/length_kkk - gama_rho24)
            rho11_new_value = rho1_1[kkk-1][i] + delt*M_rho11*f_rho_1_1
            rho21_new_value = rho2_1[kkk-1][i] + delt*M_rho21*f_rho_2_1
            rho12_new_value = rho1_2[kkk-1][i] + delt*M_rho12*f_rho_1_2
            rho22_new_value = rho2_2[kkk-1][i] + delt*M_rho22*f_rho_2_2
            rho13_new_value = rho1_3[kkk-1][i] + delt*M_rho11*f_rho_1_3
            rho23_new_value = rho2_3[kkk-1][i] + delt*M_rho21*f_rho_2_3
            rho14_new_value = rho1_4[kkk-1][i] + delt*M_rho12*f_rho_1_4
            rho24_new_value = rho2_4[kkk-1][i] + delt*M_rho22*f_rho_2_4
            
            l1_new = - (h14 * rho14_new_value + h24 * rho24_new_value) * np.cos(np.pi/2) + (h13 * rho13_new_value + h23 * rho23_new_value) * np.cos(3*np.pi/4)
            l2_new = - (h14 * rho14_new_value + h24 * rho24_new_value) * np.sin(np.pi/2) + (h13 * rho13_new_value + h23 * rho23_new_value) * np.sin(3*np.pi/4)
            f_g1 = -(M_11*tao1_k*beta13*h13*rho13_new_value + M_21*tao2_k*beta23*h23*rho23_new_value) - (stablize_l(np.sign(h13*rho13_new_value + h23*rho23_new_value)/(h14*rho14_new_value + h24*rho24_new_value))*(beta14*h14*rho14_new_value + beta24*h24*rho24_new_value)*(M_11*tao1_k1*abs(h13*rho13_new_value)+M_21*tao2_k1*abs(h23*rho23_new_value))  )
            f_g2 = -(M_12*tao1_k1*beta14*h14*rho14_new_value + M_22*tao2_k1*beta24*h24*rho24_new_value) - (stablize_l(np.sign(h14*rho14_new_value + h24*rho24_new_value)/(h13*rho13_new_value + h23*rho23_new_value))*(beta13*h13*rho13_new_value + beta23*h23*rho23_new_value)*(M_12*tao1_k*abs(h14*rho14_new_value)+M_22*tao2_k*abs(h24*rho24_new_value))  )
            f_g1_r = f_g1 * np.cos(2*np.pi/4) + f_g2 * np.cos(3*np.pi/4)
            f_g2_r = f_g1 * np.sin(2*np.pi/4) + f_g2 * np.sin(3*np.pi/4)
            l_1 = -(h14 * rho14_new_value + h24 * rho24_new_value)
            l_2 = (h13 * rho13_new_value + h23 * rho23_new_value)

        if l1_old != 0 and l2_old/l1_old >= -1 and l2_old/l1_old < 0:       
            gama_rho11 = gamma1*np.sign(h11*rho1_1[kkk-1][i] + h21*rho2_1[kkk-1][i])*sca*h11
            gama_rho21 = gamma1*np.sign(h11*rho1_1[kkk-1][i] + h21*rho2_1[kkk-1][i])*sca*h21
            gama_rho12 = 0
            gama_rho22 = 0
            gama_rho13 = 0
            gama_rho23 = 0
            gama_rho14 = gamma4*np.sign(h14*rho1_4[kkk-1][i] + h24*rho2_4[kkk-1][i])*sca*h14
            gama_rho24 = gamma4*np.sign(h14*rho1_4[kkk-1][i] + h24*rho2_4[kkk-1][i])*sca*h24
            f_rho_1_1 = stablize(G_r*(beta11*h11*(f_R_1_1 + f_th_1_1)) + LL_sca*lamda_new_value * (2*h11*(h11*rho1_1[kkk-1][i] + h21*rho2_1[kkk-1][i])-2*h11*(h14*rho1_4[kkk-1][i]+h24*rho2_4[kkk-1][i])*coeff1)/length_kkk - gama_rho11)
            f_rho_2_1 = stablize(G_r*(beta21*h21*(f_R_2_1 + f_th_2_1)) + LL_sca*lamda_new_value * (2*h21*(h11*rho1_1[kkk-1][i] + h21*rho2_1[kkk-1][i])-2*h21*(h14*rho1_4[kkk-1][i]+h24*rho2_4[kkk-1][i])*coeff1)/length_kkk - gama_rho21)
            f_rho_1_2 = 0
            f_rho_2_2 = 0
            f_rho_1_3 = 0
            f_rho_2_3 = 0
            f_rho_1_4 = stablize(G_r*(beta14*h14*(f_R_1_4 + f_th_1_4)) + LL_sca*lamda_new_value * (2*h14*(h14*rho1_4[kkk-1][i] + h24*rho2_4[kkk-1][i])-2*h14*(h11*rho1_1[kkk-1][i]+h21*rho2_1[kkk-1][i])*coeff1)/length_kkk - gama_rho14)
            f_rho_2_4 = stablize(G_r*(beta24*h24*(f_R_2_4 + f_th_2_4)) + LL_sca*lamda_new_value * (2*h24*(h14*rho1_4[kkk-1][i] + h24*rho2_4[kkk-1][i])-2*h24*(h11*rho1_1[kkk-1][i]+h21*rho2_1[kkk-1][i])*coeff1)/length_kkk -  gama_rho24)
            rho11_new_value = rho1_1[kkk-1][i] + delt*M_rho11*f_rho_1_1
            rho21_new_value = rho2_1[kkk-1][i] + delt*M_rho21*f_rho_2_1
            rho12_new_value = rho1_2[kkk-1][i] + delt*M_rho12*f_rho_1_2
            rho22_new_value = rho2_2[kkk-1][i] + delt*M_rho22*f_rho_2_2
            rho13_new_value = rho1_3[kkk-1][i] + delt*M_rho11*f_rho_1_3
            rho23_new_value = rho2_3[kkk-1][i] + delt*M_rho21*f_rho_2_3
            rho14_new_value = rho1_4[kkk-1][i] + delt*M_rho12*f_rho_1_4
            rho24_new_value = rho2_4[kkk-1][i] + delt*M_rho22*f_rho_2_4
            
            l1_new =  -(h11 * rho11_new_value + h21 * rho21_new_value) * np.cos(3*np.pi/4) + (h14 * rho14_new_value + h24 * rho24_new_value) * np.cos(0*np.pi)
            l2_new =  -(h11 * rho11_new_value + h21 * rho21_new_value) * np.sin(3*np.pi/4) + (h14 * rho14_new_value + h24 * rho24_new_value) * np.sin(0*np.pi)
            f_g1 = -(M_12*tao1_k*beta14*h14*rho14_new_value + M_22*tao2_k*beta24*h24*rho24_new_value) - (stablize_l(np.sign(h14*rho14_new_value + h24*rho24_new_value)/(h11*rho11_new_value + h21*rho21_new_value))*(beta11*h11*rho11_new_value + beta21*h21*rho21_new_value)*(M_12*tao1_k1*abs(h14*rho14_new_value)+M_22*tao2_k1*abs(h24*rho24_new_value))  )
            f_g2 = -(M_11*tao1_k1*beta11*h11*rho11_new_value + M_21*tao2_k1*beta21*h21*rho21_new_value) - (stablize_l(np.sign(h11*rho11_new_value + h21*rho21_new_value)/(h14*rho14_new_value + h24*rho24_new_value))*(beta14*h14*rho14_new_value + beta24*h24*rho24_new_value)*(M_11*tao1_k*abs(h11*rho11_new_value)+M_21*tao2_k*abs(h21*rho21_new_value))  )
            f_g1_r = f_g1 * np.cos(3*np.pi/4) + f_g2 * np.cos(0*np.pi/4)
            f_g2_r = f_g1 * np.sin(3*np.pi/4) + f_g2 * np.sin(0*np.pi/4)
            l_2 = (h14 * rho14_new_value + h24 * rho24_new_value)
            l_1 = -(h11 * rho11_new_value + h21 * rho21_new_value)



        rho1_1_new.append(rho11_new_value)
        rho2_1_new.append(rho21_new_value)
        rho1_2_new.append(rho12_new_value)
        rho2_2_new.append(rho22_new_value)
        rho1_3_new.append(rho13_new_value)
        rho2_3_new.append(rho23_new_value)
        rho1_4_new.append(rho14_new_value)
        rho2_4_new.append(rho24_new_value)
        
        
        eta = 1.3
        sinhsin = np.sinh(5*l_2/eta)
        sinhcos = np.sinh(5*l_1/eta)
        coshsin = np.cosh(5*l_2/eta)
        coshcos = np.cosh(5*l_1/eta)
        Rn_1_l2 = l_1*sinhsin/(1+coshsin)/2
        Rn_2_l2 = -(l_2*sinhsin/(1+coshsin) + 5*l_1**2/eta*(coshsin/(coshsin+1) - sinhsin**2/(coshsin+1)**2))/2
        Rn_1_l1 = -l_2*sinhcos/(1+coshcos)/2
        Rn_2_l1 = -(l_1*sinhcos/(coshcos+1) + 5*l_2**2/eta*(coshcos/(coshcos+1) - sinhcos**2/(coshcos+1)**2))  /2      
        GGG = 0.25 + gamma2*(Rn(l_2) + Rn_2_l2) + gamma1*(Rn(l_1) + Rn_2_l1) 


        x1_new_value = x1[kkk-1][i] + delt*((kap*GGG)*-1*l2_car+ -f_g1_r + K_sp*qta*l1_car)
        x2_new_value = x2[kkk-1][i] + delt*((kap*GGG)*l1_car+ -1*f_g2_r + K_sp*qta*l2_car) 

        x1_new.append(x1_new_value)
        x2_new.append(x2_new_value)

        if i == 0: 
            x1_0 = x1_new_value
            x2_0 = x2_new_value
        tolerance = 0.3
        if i != 0 and math.sqrt((x1_new_value-x1_new[i-1])**2+(x2_new_value-x2_new[i-1])**2) < tolerance:
            close_num.append(i)
            close += 1
    for iiii in range (0, len(close_num)):
        x1_new.pop(close_num[len(close_num)-iiii-1])
        x2_new.pop(close_num[len(close_num)-iiii-1])
    n = n - close
    if kkk%1000 == 0: print(n)
    x1_new.append(x1_0)
    x2_new.append(x2_0)
    x1.append(x1_new)
    x2.append(x2_new)
    lamda.append(lamda_new)
    rho1_1.append(rho1_1_new)
    rho1_2.append(rho1_2_new)
    rho2_1.append(rho2_1_new)
    rho2_2.append(rho2_2_new)
    rho1_3.append(rho1_3_new)
    rho1_4.append(rho1_4_new)
    rho2_3.append(rho2_3_new)
    rho2_4.append(rho2_4_new)    
    if (kkk*delt)%50 ==0:
        t_1.append(delt*kkk)
        sum_1.append(abs(np.trapz(x2_new, x1_new)/R0**2))
        x1_zero.append(x1_new[0])

    if (kkk*delt)%(100) == 0: 
        fig, ax = plt.subplots(figsize=(3, 3))
        ax.plot(x1[kkk], x2[kkk], linewidth = 2.0)
        ax.set_xlim(xmin=-150, xmax=150)
        ax.set_ylim(ymin=-150, ymax=150)
        ax.set_title('t = ' + str(kkk*delt))
        plt.tick_params(labelsize = 16)
        labels = ax.get_xticklabels() + ax.get_yticklabels()
        [label.set_fontname('Arial') for label in labels]
        plt.show()



font = {'family': 'Arial', 'weight': 'normal', 'size': 23}

fig, ax = plt.subplots(figsize=(3, 3), dpi = 600)
for iii in range (0,int(All/delt+1)):
    if iii%(int(All/delt/5)) == 0:
        ax.plot(x1[iii], x2[iii], linewidth = 2.0, c = 'black')
ax.set_xlim(xmin=-150, xmax=150)
ax.set_ylim(ymin=-150, ymax=150)
ax.set_xlabel(r'$\tilde{x}_1$', font)
ax.set_ylabel(r'$\tilde{x}_2$', font)
plt.tick_params(labelsize = 16)
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Arial') for label in labels]
plt.show()

font = {'family': 'Arial', 'weight': 'normal', 'size': 23}
fig, ax1 = plt.subplots(figsize=(5, 3), dpi = 600)
ax1.plot(t_1, np.array(sum_1)/3.2, c='black', linewidth = 2.0)
ax1.set_xlim(xmin=0, xmax=All)
ax1.set_ylim(ymin=0.5, ymax=1.0)
ax1.set_xlabel(r'$\tilde{t}$', font)
ax1.set_ylabel(r'$\tilde{S}$', font)
plt.tick_params(labelsize = 16)
labels = ax1.get_xticklabels() + ax1.get_yticklabels()
[label.set_fontname('Arial') for label in labels]
plt.show()



