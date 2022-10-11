##Plots results of study of interface width effect on diffusion mass change in a 2 grain systen --> Fig. 9 (a)

# import packages
import matplotlib.pyplot as plt
import matplotlib
import csv
import math
import numpy as np
from scipy import stats
from IPython.display import set_matplotlib_formats
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

##Time vs dM --> 1D model vs analytical expression
fig = plt.figure(figsize=(3,1.85),dpi=500)
ax = fig.add_subplot(111)


# user-defined parameters
workdir = 'workdir.'

csv_dirName = '2_grain/2_grain.csv' # directory containing csv file

params_data = readDatFile('dakota_tabular.dat')

params_data[0]+=  ["mass_loss","elapsed"]#['inital_total_cr',"final_total_cr","inital_metal_cr","final_metal_cr","corrosion_depth"]
#
for i in range(1,len(params_data)):
    cw_dir = workdir+str(i)+"/"
    data   = np.asarray(getRawData(cw_dir+csv_dirName,','))

    end_time = data[0,-1]
    if (end_time < 1e5):
        print("Incomplete simulation", (i ))
    inital_metal_cr = data[-1,1]
    final_metal_cr  = data[-1,-1]

    ax.plot(data[0,1:]/3600,data[-1,1:]-data[-1,1],marker='*',linewidth=0,markevery=10,ms=3)
    mass_change = final_metal_cr - inital_metal_cr
    params_data[i]+= [mass_change,end_time]

param_names = params_data[0]
sens_data   = np.asarray(params_data[1:])

X = sens_data[:,0]
Y = sens_data[:,1]

P = np.polyfit(X,Y,1)
y_fit = Y[-1]#np.polyval(P,[5e-4])
print("Ideal mass loss ",y_fit)

# t = np.linspace(0,1e6,1000)
# dC = 13.5*(2.07e-7*t)**(1/2) - 6.89e-4*t**(3/4)*(exp(-224.0/t*(1/4)) - 1.0) + 6.89e-4*t*(3/4)*exp(-4.34e+7/t)*(exp(-224.0/t**(1/4)) - 1.0) - 8.05*t**(1/4)*erfc(6590.0/t**(1/2))*(exp(-224.0/t**(1/4)) - 1.0)
#
# plt.plot(t/3600,dC,'k--')

err = (Y-y_fit)*100/y_fit
print(Y)

xmin = 0
xmax = 300
ymin = 0
ymax = 60

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

xticks=[ int(tick) for tick in np.linspace(xmin,xmax, 4)]
ax.set_xticks( xticks )
ax.set_xticklabels( xticks ,fontsize=6)

#
yticks = np.linspace(ymin,ymax, 4)
ax.set_yticks( yticks)
ax.set_yticklabels( np.around(yticks,1),fontsize=6)

ax.tick_params(axis='x',direction='in',length=5)
ax.tick_params(axis='y',direction='in',which='major',length=5)

#Set axis spine widths
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)


ax.margins(x=0, y=0)
ax.set_ylabel(r'$\Delta c_{Cr}(mol \cdot \mu$m$^3)$',fontsize=7,fontweight='bold')
ax.set_xlabel('Time (hrs)',fontsize=7,fontweight='bold')

legend_properties = {'weight':'bold','size':6}
lgd = ax.legend([r'$l=2\mu$m',r'$l=1\mu$m',r'$l=0.5\mu$m',r'$l=0.2\mu$m',r'$l=0.1\mu$m'],prop=legend_properties,framealpha=0) #,ncol=2


ax.margins(x=0, y=0)


fig.tight_layout(pad=0.3)  #Removes unnecessary padding from figure. Usually makes figure look much better
# plt.plot(X,Y,'k*')
# plt.show()
fig.savefig(PATH+'gb_width_diffusion.png',dpi=500, transparent=True)
