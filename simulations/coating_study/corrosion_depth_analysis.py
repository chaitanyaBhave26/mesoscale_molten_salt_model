# import packages
import matplotlib.pyplot as plt
import matplotlib
import csv
import math
import numpy as np
from scipy import stats
# from IPython.display import set_matplotlib_formats
from numpy import exp
from scipy.special import erfc

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

def get_plot_vars(training_data,idx):
    T = training_data[0]
    H = T/3600 #time in hours
    M = training_data[idx]
    dM = ((M[1] - M[1:])) #mass change in mg/cm2
    return(H,dM)

def readDatFile(fileName):
    rawData=[]
    with open(fileName,'r') as f:
        rows = f.readlines()
    for row in rows:
        line = row.split()
        del line[1]
        del line[0]
        del line[-1]
        for (idx,val) in enumerate(line):
            try:
                line[idx] = float(val)
            except:
                pass
        rawData+=[line]
    return rawData

plt.close('all')
PATH = '/blue/michael.tonks/chaitanya.bhave/yellowjacket/yj_first_paper/Figures/'
#set plot dimensions

plt.rcParams.update({'font.family':'Arial'})
plt.rc('font', family='sans-serif',weight='bold')


data_dakota = readDatFile('dakota_tabular.dat')
# data_dakota[0]+=  ["corr_depth"]

data_corr_depth = getRawData('corr_depth.txt',' ')

corr_depth = 150 - np.asarray(data_corr_depth[1])

param_names = data_dakota[0]
sens_data   = np.asarray(data_dakota[1:])

# #Effect of G_alloy
X = sens_data[:,0]
Y = corr_depth
print(X,Y)
fig = plt.figure(figsize=(3,1.85),dpi=500)
ax = fig.add_subplot(111)
ax.plot(X,Y,'k*')
P = np.polyfit(X,Y,1)

S_1 = P[0]*Y[0]/X[0]

X_plot = np.linspace(5,15,10)
ax.plot(X_plot,np.polyval(P,X_plot),'r--')

xmin = 4
xmax = 16
ymin = 15
ymax = 40

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

xticks=[ int(tick) for tick in np.linspace(xmin,xmax, 4)]
ax.set_xticks( xticks )
ax.set_xticklabels( xticks ,fontsize=6)

#
yticks = np.linspace(ymin,ymax, 6)
ax.set_yticks( yticks)
ax.set_yticklabels( np.around(yticks,2),fontsize=6)

ax.set_ylabel(r'Corrosion depth $\mu$m',fontsize=7,fontweight='bold')
ax.set_xlabel(r'Alloy grain size ($\mu$m)',fontsize=7,fontweight='bold')

fig.tight_layout(pad=0.3)  #Removes unnecessary padding from figure. Usually makes figure look much better
fig.savefig('G_alloy_corr_depth.png')
#


# #Effect of G_coating
X = sens_data[:,1] #
Y = corr_depth #
# print(X,Y)
fig = plt.figure(figsize=(3,1.85),dpi=500)
ax = fig.add_subplot(111)
ax.plot(X,Y,'k*')
P = np.polyfit(X,Y,1)
# print(P*Y[0]/X[0])
S_2 = P[0]*Y[0]/X[0]

X_plot = np.linspace(5,15,10)
ax.plot(X_plot,np.polyval(P,X_plot),'r--')

xmin = 4
xmax = 16
ymin = 15
ymax = 40

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

xticks=[ int(tick) for tick in np.linspace(xmin,xmax, 4)]
ax.set_xticks( xticks )
ax.set_xticklabels( xticks ,fontsize=6)

#
yticks = np.linspace(ymin,ymax, 6)
ax.set_yticks( yticks)
ax.set_yticklabels( np.around(yticks,2),fontsize=6)

ax.set_ylabel(r'Corrosion depth $\mu$m',fontsize=7,fontweight='bold')
ax.set_xlabel(r'Coating grain size ($\mu$m)',fontsize=7,fontweight='bold')
fig.tight_layout(pad=0.3)  #Removes unnecessary padding from figure. Usually makes figure look much better
fig.savefig('G_coating_corr_depth.png')


# #Effect of coating_thickness
X = sens_data[:,2] #
Y = corr_depth #
# print(X,Y)
fig = plt.figure(figsize=(3,1.85),dpi=500)
ax = fig.add_subplot(111)
ax.plot(X,Y,'k*')
P = np.polyfit(X,Y,1)
print(P*Y[0]/X[0])
S_3 = P[0]*Y[0]/X[0]

X_plot = np.linspace(30,90,10)
ax.plot(X_plot,np.polyval(P,X_plot),'r--')

xmin = 20
xmax = 100
ymin = 15
ymax = 40

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

xticks=[ int(tick) for tick in np.linspace(xmin,xmax, 5)]
ax.set_xticks( xticks )
ax.set_xticklabels( xticks ,fontsize=6)

#
yticks = np.linspace(ymin,ymax, 6)
ax.set_yticks( yticks)
ax.set_yticklabels( np.around(yticks,1),fontsize=6)

ax.set_ylabel(r'Corrosion depth $\mu$m',fontsize=7,fontweight='bold')
ax.set_xlabel(r'Coating thickness ($\mu$m)',fontsize=7,fontweight='bold')
fig.tight_layout(pad=0.3)  #Removes unnecessary padding from figure. Usually makes figure look much better
fig.savefig('coating_thickness_corr_depth.png')

X=[1,2,3] #,3
S = [S_1,S_2,S_3] #,S_3
S = S/max(np.abs(S))
print("Sensitivity (normalized) = ",S)

fig = plt.figure(figsize=(3,1.85),dpi=500)
ax = fig.add_subplot(111)
ax.bar(X,S,color='b',alpha=0.8)
ax.axhline(y=0,linestyle='--',linewidth=1,color='k')

ax.set_xticks(X)
ax.set_xticklabels( [r"$G_{alloy}$",r"$G_{coating}$",r"$T_{coating}$"],fontsize=6) #,r"$T_{coating}$"

yticks = np.linspace(-1,1, 5)
ax.set_yticks( yticks)
ax.set_yticklabels( np.around(yticks,1),fontsize=6)

ax.set_ylabel('Normalized sensitivities',fontsize=7,fontweight='bold')
ax.set_xlabel('Sampling parameter',fontsize=7,fontweight='bold')

fig.tight_layout(pad=0.3)  #Removes unnecessary padding from figure. Usually makes figure look much better
fig.savefig("Sensitivity_bar_plot_Corr_depth.png")