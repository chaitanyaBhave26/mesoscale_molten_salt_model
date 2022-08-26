##Used to fit an ideal solution free energy to the Ni-Cr CALPHAD assessment data

# import packages
import matplotlib.pyplot as plt
import csv
import math
import numpy as np
from scipy import stats
import scipy.optimize as optimization

PATH = '/blue/michael.tonks/chaitanya.bhave/yellowjacket/yj_first_paper/Figures/'

def getRawData(fileName, delim): # extract raw data from csv file

    rawData = []
    with open(fileName, 'r') as f:
        CSVReader = csv.reader(f, delimiter = delim, skipinitialspace = True)
        array_size = len(next(CSVReader))
        rawData = [ [] for i in range(array_size)]
        for row in CSVReader:
            for (i,val) in enumerate(row):
                rawData[i].append(float(val) )
    return rawData

def idealFreeEnergy(xdata,H_Cr_metal,S_Cr_metal):
    x_Ni = xdata[0,:]
    T = xdata[1,:]
    G0_f_Ni = -5179.159 + 117.854*T - 22.096*T*np.log(T) - (4.8407e-3)*T**2 #J/mol
    G0_f_Cr = -1572.94 + 157.643*T - 26.908*T*np.log(T) + 1.89435E-3*T**2 - 1.47721E-6*T**3 + 139250/T #J/mol
    R = 8.3145
    F = G0_f_Ni*x_Ni + (G0_f_Cr + H_Cr_metal - S_Cr_metal*T) *(1-x_Ni) + R*T*(x_Ni*np.log(x_Ni) + (1-x_Ni)*np.log(1-x_Ni))
    return F/1e3


def regularFreeEnergy(x_Ni,E0_Cr_metal,beta):
    return (E0_Cr_metal*(1-x_Ni) + 8.314*1000*(x_Ni*np.log(x_Ni + 1e-12) + (1-x_Ni)*np.log(1-x_Ni + 1e-12)  ) + beta*x_Ni*(1-x_Ni) )

def CalphadEnergy(T,x_Ni):
    G0_f_Ni = -5179.159 + 117.854*T - 22.096*T*np.log(T) - (4.8407e-3)*T**2 #J/mol
    G0_f_Cr = -1572.94 + 157.643*T - 26.908*T*np.log(T) + 1.89435E-3*T**2 - 1.47721E-6*T**3 + 139250/T #J/mol
    R = 8.3145 #J/mol/K

    T_C = -1109*(1-x_Ni) + 633*x_Ni -3605*x_Ni*(1-x_Ni)
    T_C = T_C*(T_C>0) +  T_C*(T_C<0)/-3

    B0 = -2.46*(1-x_Ni) +0.528*x_Ni - 1.91*(1-x_Ni)*x_Ni
    B0 = B0*(B0>0) +  B0*(B0<0)/-3

    p = 0.28
    tau = T/T_C
    A = (518/1125) + (11692/15975)*( (1/p) - 1)
    f_tau = (1/(10*tau**5) + 1/(315*tau**15) + 1/(1500*tau**25)  )/A #Our system is always above T_C

    G_mag = R*T*np.log(1 + B0)*f_tau

    L = 8347 -12.1038*T + (29895 - 16.3838*T)*(1-2*x_Ni)

    F = G0_f_Ni*x_Ni + G0_f_Cr*(1-x_Ni) + R*T*(x_Ni*np.log(x_Ni) + (1-x_Ni)*np.log(1-x_Ni) ) + L*x_Ni*(1-x_Ni)  + G_mag
    return F/1e3


#SGTE formation energies
T = np.array([700,800,900,1000]) #K
P = 101325 #Pa

x_Ni = np.arange(0.9,1-1e-9,1e-3)
size = (x_Ni.shape)[0]
print("Size of x_Ni is = ",size)

F_0 = CalphadEnergy(T[0],x_Ni)
F_1 = CalphadEnergy(T[1],x_Ni)
F_2 = CalphadEnergy(T[2],x_Ni)
F_3 = CalphadEnergy(T[3],x_Ni)


#Stack data points
F = np.hstack([F_0,F_1,F_2,F_3])
X = np.hstack([x_Ni,x_Ni,x_Ni,x_Ni])
T_s = np.hstack([np.repeat(temp,x_Ni.shape[0]) for temp in T ] )

xdata = np.vstack([X,T_s])

#Perform Least squares fit
P0 = np.array([0.0,0.0])
P0,cov = optimization.curve_fit(idealFreeEnergy,xdata,F,P0)
print(P0)

F_fit = idealFreeEnergy(xdata,P0[0],P0[1])

F0_fit = F_fit[:size]
F1_fit = F_fit[size:2*size]
F2_fit = F_fit[2*size:3*size]
F3_fit = F_fit[3*size:]

rmse_f0 = math.sqrt(np.sum(np.power(F_0 - F0_fit,2))/F0_fit.shape[0])
print(rmse_f0/np.mean(F_0))

# close existing plots
plt.close('all')
#set plot dimensions
fig = plt.figure(figsize=(3.0,1.85),dpi=500)
plt.rcParams.update({'font.family':'Arial'})
plt.rc('font', family='sans-serif',weight='bold')

ax = fig.add_subplot(1,1,1)

plt.xlim(0.9,1.0)
plt.ylim(-50.0,-20.0)
plt.xticks(np.arange(0.9,1.02, 0.02),fontsize=6)  # arbitrary chosen
plt.yticks(np.linspace(-50.0,-20.0, 4),fontsize=6)  # arbitrary chosen

ax.tick_params(axis='x',direction='in')#,pad=6,)
ax.tick_params(axis='y',direction='in')#,pad=6)

n = 10
plt.plot(x_Ni[1::n],F_0[1::n],'r*',markersize=4)
plt.plot(x_Ni[1::n],F_1[1::n],'b*',markersize=4)
plt.plot(x_Ni[1::n],F_2[1::n],'g*',markersize=4)
plt.plot(x_Ni[1::n],F_3[1::n],'k*',markersize=4)

plt.plot(x_Ni,F0_fit,'r--',linewidth=1.0)
plt.plot(x_Ni,F1_fit,'b--',linewidth=1.0)
plt.plot(x_Ni,F2_fit,'g--',linewidth=1.0)
plt.plot(x_Ni,F3_fit,'k--',linewidth=1.0)

legend_properties = {'weight':'bold','size':6}
lgd = plt.legend(["T = 700 K","T = 800 K","T = 900 K","T = 1000 K"], frameon=False, ncol=2,prop=legend_properties)

plt.margins(x=0, y=0)
plt.xlabel("Atomic fraction of Ni (Ni-Cr)",fontsize=6,fontweight='bold')#,labelpad=15)
plt.ylabel(r"Gibbs energy ($\times10^3$ J/mol)",fontsize=6,fontweight='bold')#,labelpad=15)

#Set axis spine widths
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)

fig.tight_layout()
fig.savefig(PATH+'Cr_excess_energy.png',dpi=500, transparent=True,bbox_inches='tight')
