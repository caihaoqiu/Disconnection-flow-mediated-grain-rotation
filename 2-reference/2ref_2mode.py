
import numpy as np
import matplotlib.pyplot as plt

def mo(x1, x2):
    length = np.sqrt(x1**2+x2**2)
    return length

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

def abs_v(list_1):
    l = []
    for i in list_1:
        l.append(abs(i))
    return l
def stablize(x):
    ee = 0.03
    return ee*np.tanh(x/ee)

e1 = [1,0]
e2 = [0,1]




n = 60
theta_step = np.pi/180
x1_0 = []
x2_0 = []
x1 = []
x2 = []
theta = []
t = [0]
R = 100
ep = 0.1
k_sp = 0.01

M_1 = 1
M_2 = 1
Ml = -1000
M_rho1 = 0.01
M_rho2 = 0.01


gamma1 = 1
gamma2 = 2


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



Ra = abs(4*h11*beta11)


G = 1e-5
a = 0.01

for ii in range(0, n):
    theta.append(360/n*ii*theta_step)
    x1_0.append(R*np.cos(theta[ii]))
    x2_0.append(R*np.sin(theta[ii]))
x1_0.append(100)
x2_0.append(0)
x1.append(x1_0)
x2.append(x2_0)
kap = []
l1_0 = []
l2_0 = []
rho1_1_0 = []
rho1_2_0 = []
rho2_1_0 = []
rho2_2_0 = []
length_0 = []
lamda_0 = []
ra = 1/3
sca = 1e-5
LL = 1e-3
for i in range(0, n):
    dev1_1 = dev(1, i, n)[1]
    dev1_2 = dev(1, i, n)[2]
    dev2_1 = dev(1, i, n)[3]
    dev2_2 = dev(1, i, n)[4]    
    ll1 = dev1_1/np.sqrt(dev1_1**2+dev1_2**2)
    ll2 = dev1_2/np.sqrt(dev1_1**2+dev1_2**2)
    rho11 = ra*ll2 / h11
    rho12 = -ra*ll1 / h12
    rho21 = (1-ra)*ll2 / h21
    rho22 = -(1-ra)*ll1 / h22
    l1_0.append(ll1)
    l2_0.append(ll2)       
    rho1_1_0.append(rho11)        #\rho_1^{(1)}
    rho1_2_0.append(rho12)        #\rho_1^{(2)}    
    rho2_1_0.append(rho21)        #\rho_2^{(1)}
    rho2_2_0.append(rho22)        #\rho_2^{(2)}    
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
lamda = [lamda_0]

rho1_1_v = [sum(abs_v(rho1_1_0))]
rho1_2_v = [sum(abs_v(rho1_2_0))]
rho2_1_v = [sum(abs_v(rho2_1_0))]
rho2_2_v = [sum(abs_v(rho2_2_0))]
dxdx =[]
for ll in range(0, len(x1_0)-1):
    if ll == 0:
        dxdx.append(np.sqrt((x1_0[len(x1_0)-2] - x1_0[1])**2 + (x2_0[len(x1_0)-2] - x2_0[1])**2)/2)
    if ll == len(x1_0)-2:
        dxdx.append(np.sqrt((x1_0[ll-1] - x1_0[0])**2 + (x2_0[ll-1] - x2_0[0])**2)/2)
    if ll < len(x1_0)-2 and ll > 0:
        dxdx.append(np.sqrt((x1_0[ll-1] - x1_0[ll+1])**2 + (x2_0[ll-1] - x2_0[ll+1])**2)/2)       
rho1_1_vv = [sum(np.multiply(np.array(abs_v(rho1_1_0)), np.array(dxdx)))/sum(dxdx)]
rho1_2_vv = [sum(np.multiply(np.array(abs_v(rho1_2_0)), np.array(dxdx)))/sum(dxdx)]
rho2_1_vv = [sum(np.multiply(np.array(abs_v(rho2_1_0)), np.array(dxdx)))/sum(dxdx)]
rho2_2_vv = [sum(np.multiply(np.array(abs_v(rho2_2_0)), np.array(dxdx)))/sum(dxdx)]


t_1 = [0]
sum_1 = []
sum_1.append(abs(np.trapz(x2_0, x1_0)/R**2))

All = 4000
delt = 0.001


for kkk in range(1,int(All/delt+1)):
    close = 0
    close_num = []


    t.append(delt*kkk)
    x1_new = []
    x2_new = []
    lamda_new = []
    rho1_1_new = []
    rho1_2_new = []
    rho2_1_new = []
    rho2_2_new = []
    length_new = []
    kappa1 = []
    for i in range(0, n):
        qta = dev(kkk, i, n)[0]
        dev1_1kkk = dev(kkk, i, n)[1]
        dev1_2kkk = dev(kkk, i, n)[2]
        length_kkk = np.sqrt(dev1_1kkk**2 + dev1_2kkk**2)
        l1_car_old = -h12 * rho1_2[kkk-1][i] - h22 * rho2_2[kkk-1][i]
        l2_car_old =  h11 * rho1_1[kkk-1][i] + h21 * rho2_1[kkk-1][i]
        lamda_new_value = lamda[kkk-1][i] + Ml*delt/length_kkk * (l1_car_old**2 + l2_car_old**2 - 1)
        lamda_new.append(lamda_new_value)

        I12_11_all = 0
        I12_12_all = 0
        I12_21_all = 0
        I12_22_all = 0
        f_R_1_1 = 0
        f_R_2_1 = 0
        f_R_1_2 = 0
        f_R_2_2 = 0
        f_th_1_1 = 0
        f_th_2_1 = 0
        f_th_1_2 = 0
        f_th_2_2 = 0

        for j in range(0,n):           
            dx1 = x1[kkk-1][i]-x1[kkk-1][j]
            dx2 = x2[kkk-1][i]-x2[kkk-1][j]
            RR = np.sqrt(dx1**2 + dx2**2)
            if i != j:
                costheta1 = dx2/RR
                costheta2 = -dx1/RR
            dev1_1jjj = devj(kkk, j, n)[1]
            dev1_2jjj = devj(kkk, j, n)[2]
            length_jjj = np.sqrt(dev1_1jjj**2 + dev1_2jjj**2)

            q1 = (dx1)**2 + (dx2)**2 + a**2          
            tao12_1 = (dx1/q1)*(1-2*(dx2)**2/q1)
            tao12_2 = (dx2/q1)*(1-2*(dx2)**2/q1)

            I12_11_all += h11*rho1_1[kkk-1][j]*length_jjj * tao12_1
            I12_12_all += h12*rho1_2[kkk-1][j]*length_jjj * tao12_2
            I12_21_all += h21*rho2_1[kkk-1][j]*length_jjj * tao12_1
            I12_22_all += h22*rho2_2[kkk-1][j]*length_jjj * tao12_2
            if i != j:
                f_R_1_1 += (2*rho1_1[kkk-1][j]*beta11*h11 + rho2_1[kkk-1][j]*beta21*h21)*np.log(RR/Ra)*length_jjj
                f_R_2_1 += (rho1_1[kkk-1][j]*beta11*h11 + 2*rho2_1[kkk-1][j]*beta21*h21)*np.log(RR/Ra)*length_jjj
                f_R_1_2 += (2*rho1_2[kkk-1][j]*beta12*h12 + rho2_2[kkk-1][j]*beta22*h22)*np.log(RR/Ra)*length_jjj
                f_R_2_2 += (rho1_2[kkk-1][j]*beta12*h12 + 2*rho2_2[kkk-1][j]*beta22*h22)*np.log(RR/Ra)*length_jjj

                f_th_1_1 += (2*rho1_1[kkk-1][j]*beta11*h11*costheta1*costheta1 + rho2_1[kkk-1][j]*beta21*h21*costheta1*costheta1 + rho1_2[kkk-1][j]*beta12*h12*costheta1*costheta2 + rho2_2[kkk-1][j]*beta22*h22*costheta1*costheta2)*length_jjj
                f_th_2_1 += (rho1_1[kkk-1][j]*beta11*h11*costheta1*costheta1 + 2*rho2_1[kkk-1][j]*beta21*h21*costheta1*costheta1 + rho1_2[kkk-1][j]*beta12*h12*costheta1*costheta2 + rho2_2[kkk-1][j]*beta22*h22*costheta1*costheta2)*length_jjj
                f_th_1_2 += (rho1_1[kkk-1][j]*beta11*h11*costheta1*costheta2 + rho2_1[kkk-1][j]*beta21*h21*costheta1*costheta2 + 2*rho1_2[kkk-1][j]*beta12*h12*costheta2*costheta2 + rho2_2[kkk-1][j]*beta22*h22*costheta2*costheta2)*length_jjj
                f_th_2_2 += (rho1_1[kkk-1][j]*beta11*h11*costheta1*costheta2 + rho2_1[kkk-1][j]*beta21*h21*costheta1*costheta2 + rho1_2[kkk-1][j]*beta12*h12*costheta2*costheta2 + 2*rho2_2[kkk-1][j]*beta22*h22*costheta2*costheta2)*length_jjj              
        tao_all = G*(beta11*I12_11_all + beta12*I12_12_all + beta21*I12_21_all + beta22*I12_22_all)

        gama_rho11 = gamma1*np.sign(l2_car_old)*sca*h11
        gama_rho12 = gamma2*np.sign(-l1_car_old)*sca*h12
        gama_rho21 = gamma1*np.sign(l2_car_old)*sca*h21
        gama_rho22 = gamma2*np.sign(-l1_car_old)*sca*h22

        
        f_rho_1_1 = stablize((G*(beta11*h11*(f_R_1_1 + f_th_1_1)) + (2*LL*h11*lamda_new_value)/length_kkk * (l2_car_old)) - gama_rho11)
        f_rho_2_1 = stablize((G*(beta21*h21*(f_R_2_1 + f_th_2_1)) + (2*LL*h21*lamda_new_value)/length_kkk * (l2_car_old)) - gama_rho21)
        f_rho_1_2 = stablize((G*(beta12*h12*(f_R_1_2 + f_th_1_2)) + (2*LL*h12*lamda_new_value)/length_kkk * (-l1_car_old)) - gama_rho12)
        f_rho_2_2 = stablize((G*(beta22*h22*(f_R_2_2 + f_th_2_2)) + (2*LL*h22*lamda_new_value)/length_kkk * (-l1_car_old)) - gama_rho22)


        rho11_new_value = rho1_1[kkk-1][i] + delt*M_rho1*f_rho_1_1
        rho21_new_value = rho2_1[kkk-1][i] + delt*M_rho2*f_rho_2_1
        rho12_new_value = rho1_2[kkk-1][i] + delt*M_rho1*f_rho_1_2
        rho22_new_value = rho2_2[kkk-1][i] + delt*M_rho2*f_rho_2_2


        rho1_1_new.append(rho11_new_value)
        rho2_1_new.append(rho21_new_value)
        rho1_2_new.append(rho12_new_value)
        rho2_2_new.append(rho22_new_value)
        
        length_new_value = length[kkk-1][i] + delt * gamma1/length[kkk-1][i]
        length_new.append(length_new_value)
        
        l1_car = -h12 * rho12_new_value - h22 * rho22_new_value 
        l2_car =  h11 * rho11_new_value + h21 * rho21_new_value 
        
        Gamma = 0.5 + 2*ep/np.pi*(gamma2*l2_car**2/(ep**2+l1_car**2)+gamma1*l1_car**2/(ep**2+l2_car**2))

        l1_1 = dev(kkk, i, n)[1]
        l2_1 = dev(kkk, i, n)[2]
        dl1_car = dev(kkk, i, n)[3]
        dl2_car = dev(kkk, i, n)[4]  
        l1_carl = l1_1/np.sqrt(l1_1**2 + l2_1**2)
        l2_carl = l2_1/np.sqrt(l1_1**2 + l2_1**2)
        kappa = (l1_1*dl2_car-dl1_car*l2_1)/(l1_1**2+l2_1**2)**1.5
        
        f_g1 = tao_all*(beta11*h11*rho1_1[kkk-1][i] + (beta12*h12*rho1_2[kkk-1][i]+beta22*h22*rho2_2[kkk-1][i])*np.tanh(np.sign(l2_car)/l1_car)*abs(h11*rho1_1[kkk-1][i]) + beta21*h21*rho2_1[kkk-1][i] + (beta12*h12*rho1_2[kkk-1][i]+beta22*h22*rho2_2[kkk-1][i])*np.tanh(np.sign(l2_car)/l1_car)*abs(h21*rho2_1[kkk-1][i]))
        f_g2 = -tao_all*(beta12*h12*rho1_2[kkk-1][i] + (beta11*h11*rho1_1[kkk-1][i]+beta21*h21*rho2_1[kkk-1][i])*np.tanh(np.sign(l1_car)/l2_car)*abs(h12*rho1_2[kkk-1][i]) + beta22*h22*rho2_2[kkk-1][i] + (beta11*h11*rho1_1[kkk-1][i]+beta21*h21*rho2_1[kkk-1][i])*np.tanh(np.sign(l1_car)/l2_car)*abs(h22*rho2_2[kkk-1][i]))

        x1_new_value = x1[kkk-1][i] +  delt*(f_g1+ (Gamma*kappa)*-l2_carl) + delt*k_sp*qta*l1_carl
        x2_new_value = x2[kkk-1][i] +  delt*(f_g2+ (Gamma*kappa)*l1_carl) + delt*k_sp*qta*l2_carl
        
        x1_new.append(x1_new_value)
        x2_new.append(x2_new_value)

        if i == 0: 
            x11 = x1_new_value
            x22 = x2_new_value    
        tolerance = 0.3
        if i != 0 and (np.sqrt((x1_new_value-x1_new[i-1])**2+(x2_new_value-x2_new[i-1])**2) < tolerance):
            close_num.append(i)
            close += 1
    if kkk%(int(All/delt/10)) ==0:
        t_1.append(delt*kkk)
        sum_1.append(abs(np.trapz(x2_new, x1_new)/R**2))   
        rho1_1_v.append(sum(abs_v(rho1_1_new)))
        rho1_2_v.append(sum(abs_v(rho1_2_new)))
        rho2_1_v.append(sum(abs_v(rho2_1_new)))
        rho2_2_v.append(sum(abs_v(rho2_2_new)))
        dxdx =[]
        for ll in range(0, len(x1_new)):
                if ll == 0:
                    dxdx.append(np.sqrt((x1_new[len(x1_new)-1] - x1_new[1])**2 + (x2_new[len(x1_new)-1] - x2_new[1])**2)/2)
                if ll == len(x1_new)-1:
                    dxdx.append(np.sqrt((x1_new[ll-1] - x1_new[0])**2 + (x2_new[ll-1] - x2_new[0])**2)/2)
                if ll != len(x1_new)-1 and ll != 0:
                    dxdx.append(np.sqrt((x1_new[ll-1] - x1_new[ll+1])**2 + (x2_new[ll-1] - x2_new[ll+1])**2)/2)       
        rho1_1_vv.append(sum(np.multiply(np.array(abs_v(rho1_1_new)), np.array(dxdx)))/sum(dxdx))
        rho1_2_vv.append(sum(np.multiply(np.array(abs_v(rho1_2_new)), np.array(dxdx)))/sum(dxdx))
        rho2_1_vv.append(sum(np.multiply(np.array(abs_v(rho2_1_new)), np.array(dxdx)))/sum(dxdx))
        rho2_2_vv.append(sum(np.multiply(np.array(abs_v(rho2_2_new)), np.array(dxdx)))/sum(dxdx))


            
    for iiii in range (0, len(close_num)):
       x1_new = list(x1_new)
       x2_new = list(x2_new)
       pp = close_num[len(close_num)-iiii-1]
       
       if pp != 0:
           rho1_1_new[pp-1] = rho1_1_new[pp-1] + rho1_1_new[pp]
           rho2_1_new[pp-1] = rho2_1_new[pp-1] + rho2_1_new[pp]
           rho1_2_new[pp-1] = rho1_2_new[pp-1] + rho1_2_new[pp]
           rho2_2_new[pp-1] = rho2_2_new[pp-1] + rho2_2_new[pp]
       x1_new.pop(pp)
       x2_new.pop(pp)
       rho1_1_new.pop(pp)
       rho2_1_new.pop(pp)
       rho1_2_new.pop(pp)
       rho2_2_new.pop(pp)
       length_new.pop(pp)
       lamda_new.pop(pp)
       n = n - 1       
       
    x1_new.append(x11)    
    x2_new.append(x22)           
    lamda.append(lamda_new)
    rho1_1.append(rho1_1_new)
    rho1_2.append(rho1_2_new)
    rho2_1.append(rho2_1_new)
    rho2_2.append(rho2_2_new)
    
    length.append(length_new)
    x1.append(x1_new)
    x2.append(x2_new)
    font = {'family': 'Arial', 'weight': 'normal', 'size': 23}
    if (kkk*delt)%100 == 0: 
        fig, ax = plt.subplots(figsize=(3, 3))
        ax.plot(x1[kkk], x2[kkk], linewidth = 2.0)
        ax.set_xlim(xmin=-150, xmax=150)
        ax.set_ylim(ymin=-150, ymax=150)
        ax.set_xlabel(r'$\tilde{x}_1$', font)
        ax.set_ylabel(r'$\tilde{x}_2$', font)
        ax.set_title('t = ' + str(kkk*delt))
        plt.tick_params(labelsize = 16)
        labels = ax.get_xticklabels() + ax.get_yticklabels()
        [label.set_fontname('Arial') for label in labels]
        plt.show()


font = {'family': 'Arial', 'weight': 'normal', 'size': 23}
fig, ax = plt.subplots(figsize=(6, 6))
for iii in range (0,int(All/delt+1)):
    if iii%(int(All/delt/10)) == 0:
        ax.plot(x1[iii], x2[iii])
ax.set_xlim(xmin=-110, xmax=110)
ax.set_ylim(ymin=-110, ymax=110)
plt.tick_params(labelsize = 20)
labels = ax.get_xticklabels() + ax.get_yticklabels()
[label.set_fontname('Arial') for label in labels]
plt.show()

font = {'family': 'Arial', 'weight': 'normal', 'size': 23}
fig, ax1 = plt.subplots(figsize=(5, 3), dpi = 600)
ax1.plot(t_1, np.array(sum_1)/3.2, c='black', linewidth = 2.0)
ax1.set_xlim(xmin=0, xmax= All)
ax1.set_ylim(ymin=0, ymax=1)
ax1.set_xticks([0,All/4, All/2, All*3/4, All])
ax1.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
ax1.set_xlabel(r'$\tilde{t}$', font)
ax1.set_ylabel(r'$\tilde{A}$', font)
plt.tick_params(labelsize = 16)
labels = ax1.get_xticklabels() + ax1.get_yticklabels()
[label.set_fontname('Arial') for label in labels]
plt.show()


