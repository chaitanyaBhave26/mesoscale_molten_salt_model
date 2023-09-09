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

# user-defined parameters
workdir = 'workdir.'

csv_dirName = 'ni20_corr/ni20_corr.csv' # directory containing csv file

params_data = readDatFile('dakota_tabular.dat')

params_data[0]+=  ["mass_loss","elapsed"]#['inital_total_cr',"final_total_cr","inital_metal_cr","final_metal_cr","corrosion_depth"]
#
for i in range(1,len(params_data)):
    cw_dir = workdir+str(i)+"/"
    data   = np.asarray(getRawData(cw_dir+csv_dirName,','))

    end_time = data[0,-1]
    if (end_time < 1e5):
        print("Incomplete simulation", (i ))
    inital_metal_cr = data[2,1]
    final_metal_cr  = data[2,-1]

    # ax.plot(data[0,1:]/3600,data[-1,1:]-data[-1,1],marker='*',linewidth=0,markevery=10,ms=3)
    mass_change = inital_metal_cr - final_metal_cr
    params_data[i]+= [mass_change,end_time]

param_names = params_data[0]
sens_data   = np.asarray(params_data[1:])

#Slices are 30-50,50-70,70-90. Include floor exclude ceiling
slice1 = np.where( np.logical_and(sens_data[:,2]>=30, sens_data[:,2]<50) )[0]
slice2 = np.where( np.logical_and(sens_data[:,2]>=50, sens_data[:,2]<70) )[0]
slice3 = np.where( np.logical_and(sens_data[:,2]>=70, sens_data[:,2]<90) )[0]


# #Effect of G_alloy
X = sens_data[:,0]
Y = sens_data[:,-2]/30
# print(X,Y)
fig = plt.figure(figsize=(3,1.85),dpi=500)
ax = fig.add_subplot(111)

S_1 = []

ax.plot(X[slice1],Y[slice1],'r*',markersize=3)
P = np.polyfit(X[slice1],Y[slice1],1)
X_plot = np.linspace(5,15,10)
ax.plot(X_plot,np.polyval(P,X_plot),'r--',linewidth=1.5)
S_1 += [P[0]]

ax.plot(X[slice2],Y[slice2],'b*',markersize=3)
P = np.polyfit(X[slice2],Y[slice2],1)
X_plot = np.linspace(5,15,10)
ax.plot(X_plot,np.polyval(P,X_plot),'b--',linewidth=1.5)
S_1 += [P[0]]

ax.plot(X[slice3],Y[slice3],'g*',markersize=3)
P = np.polyfit(X[slice3],Y[slice3],1)
X_plot = np.linspace(5,15,10)
ax.plot(X_plot,np.polyval(P,X_plot),'g--',linewidth=1.5)
S_1 += [P[0]]

S_1 = np.asarray(S_1)*np.mean(Y)/np.mean(X)
print(S_1)

# P = np.polyfit(X,Y,1)
# print(P*Y[0]/X[0])
# S_1 = P[0]*Y[0]/X[0]

# X_plot = np.linspace(5,15,10)
# ax.plot(X_plot,np.polyval(P,X_plot),'r--')

xmin = 4
xmax = 16
ymin = -0.1
ymax = 0.5

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

xticks=[ int(tick) for tick in np.linspace(xmin,xmax, 4)]
ax.set_xticks( xticks )
ax.set_xticklabels( xticks ,fontsize=6)

# #
yticks = np.linspace(ymin,ymax, 7)
ax.set_yticks( yticks)
ax.set_yticklabels( np.around(yticks,2),fontsize=6)

ax.set_ylabel('Mass loss $(mg/cm^2)$',fontsize=7,fontweight='bold')

# ax.set_ylabel(r'$dM_{Cr} (mg/cm^2)$',fontsize=7,fontweight='bold')
ax.set_xlabel(r'Alloy grain size ($\mu$m)',fontsize=7,fontweight='bold')

#Set axis spine widths
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)

fig.tight_layout(pad=0.3)  #Removes unnecessary padding from figure. Usually makes figure look much better
fig.savefig('G_alloy_mass_loss.png')


#
# #Effect of G_coating
X = sens_data[:,1]
Y = sens_data[:,-2]/30 #
# print(X,Y)
fig = plt.figure(figsize=(3,1.85),dpi=500)
ax = fig.add_subplot(111)
S_2 = []

ax.plot(X[slice1],Y[slice1],'r*',markersize=3)
P = np.polyfit(X[slice1],Y[slice1],1)
X_plot = np.linspace(5,15,10)
ax.plot(X_plot,np.polyval(P,X_plot),'r--',linewidth=1.5)
S_2 += [P[0]]

ax.plot(X[slice2],Y[slice2],'b*',markersize=3)
P = np.polyfit(X[slice2],Y[slice2],1)
X_plot = np.linspace(5,15,10)
ax.plot(X_plot,np.polyval(P,X_plot),'b--',linewidth=1.5)
S_2 += [P[0]]

ax.plot(X[slice3],Y[slice3],'g*',markersize=3)
P = np.polyfit(X[slice3],Y[slice3],1)
X_plot = np.linspace(5,15,10)
ax.plot(X_plot,np.polyval(P,X_plot),'g--',linewidth=1.5)
S_2 += [P[0]]

S_2 = np.asarray(S_2)*np.mean(Y)/np.mean(X)
print(S_2)



xmin = 4
xmax = 16
ymin = -0.1
ymax = 0.5

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

xticks=[ int(tick) for tick in np.linspace(xmin,xmax, 4)]
ax.set_xticks( xticks )
ax.set_xticklabels( xticks ,fontsize=6)

# #
yticks = np.linspace(ymin,ymax, 7)
ax.set_yticks( yticks)
ax.set_yticklabels( np.around(yticks,2),fontsize=6)

#Set axis spine widths
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)

ax.set_ylabel('Mass loss $(mg/cm^2)$',fontsize=7,fontweight='bold')
ax.set_xlabel(r'Coating grain size ($\mu$m)',fontsize=7,fontweight='bold')
fig.tight_layout(pad=0.3)  #Removes unnecessary padding from figure. Usually makes figure look much better
fig.savefig('G_coating_mass_loss.png')

# #Effect of coating_thickness
X = sens_data[:,2]
Y = sens_data[:,-2]/30 #
# print(X,Y)
fig = plt.figure(figsize=(3,1.85),dpi=500)
ax = fig.add_subplot(111)
ax.plot(X,Y,'k*',markersize=3)

# min_y = min(Y)
# Y = Y - min_y + 1e-10
# P = np.polyfit(np.exp(X),Y ,1)
# print("Logscale fit",P)
# # print(P*Y[0]/X[0])
# # S_3 = P[0]*Y[0]/X[0]

# # ax.set_yscale('log')
X_plot = np.linspace(30,90,10)
# ax.plot(X_plot,np.polyval(P,np.exp(X_plot)),'b--')

m0 = 2*8.57e3*0.2*math.sqrt(0.5e-4*3.6e6*1e-8/math.pi)
print("m0 = ",m0)
Y_fit = 3.2*np.exp(-(X_plot/math.sqrt(0.5e-4*3.6e6)) ) #0.016293911836791
# Y_fit = 3.2*np.exp(-(X_plot/math.sqrt(0.016293911836791*3.6e6)) ) #0.016293911836791
ax.plot(X_plot,Y_fit,'r--')

xmin = 20
xmax = 100
ymin = -0.1
ymax = 0.5

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

xticks=[ int(tick) for tick in np.linspace(xmin,xmax, 5)]
ax.set_xticks( xticks )
ax.set_xticklabels( xticks ,fontsize=6)

#
yticks = np.linspace(ymin,ymax, 7)
ax.set_yticks( yticks)
ax.set_yticklabels( np.around(yticks,1),fontsize=6)

ax.set_ylabel(r'$dM_{Cr} (mg/cm^2)$',fontsize=7,fontweight='bold')
ax.set_xlabel(r'Coating thickness ($\mu$m)',fontsize=7,fontweight='bold')
fig.tight_layout(pad=0.3)  #Removes unnecessary padding from figure. Usually makes figure look much better
fig.savefig('coating_thickness_mass_loss.png')

S = np.hstack([S_1,S_2])
S = S/max(np.abs(S))
print(S)

fig = plt.figure(figsize=(3,1.85),dpi=500)
ax = fig.add_subplot(111)
plt.rcParams['legend.title_fontsize'] = 6

X = ["Alloy grain size","Coating grain size"]
x_axis = np.arange(len(X))

ax.bar(x_axis-0.25,[S[0],S[3]],color='r',width = 0.2, label = "30-50")
ax.bar(x_axis,[S[1],S[4]],color='b',width = 0.2, label = "50-70")
ax.bar(x_axis+0.25,[S[2],S[5]],color='g',width = 0.2, label = "70-90")

# ax.bar([1,2,3,4,5,6],S,color=['r','g','b','r','g','b'])
# ax.bar([1,1,1,2,2,2],S,color=['r','g','b','r','g','b'])

ymin = -1.1
ymax = 0.1

# ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

xticks=[0,1]
ax.set_xticks( xticks )
ax.set_xticklabels( ["$G_{alloy}$","$G_{coating}$"] ,fontsize=6)

# #
yticks = np.linspace(0.0,-1, 5)
ax.set_yticks( yticks)
ax.set_yticklabels( np.around(yticks,2),fontsize=6)



ax.set_xlabel('Grain size parameters',fontsize=7,fontweight='bold')
ax.set_ylabel('Normalized sensitivities',fontsize=7,fontweight='bold')

legend_properties = {'weight':'bold','size':6}
lgd = ax.legend(["30-50","50-70","70-90"],prop=legend_properties,framealpha=0,title='Coating thickness ($\mu$m)')

# fig.legend(["1","2","3"])  #(["-3SD -- -SD","-SD -- SD","SD -- 3SD"])

fig.tight_layout(pad=0.3)  #Removes unnecessary padding from figure. Usually makes figure look much better
fig.savefig('Grain_size_sensitivity_mass_loss.png')



# S = [S_1,S_2,S_3]
# S = S/max(np.abs(S))
# print("Sensitivity (normalized) = ",S)
