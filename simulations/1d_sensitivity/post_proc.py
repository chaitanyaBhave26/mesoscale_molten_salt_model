# import packages
import matplotlib.pyplot as plt
import matplotlib
import csv
import math
import numpy as np
from scipy import stats
from IPython.display import set_matplotlib_formats

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
    # return(H,M[1:])
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


# user-defined parameters
workdir = 'workdir.'

csv_dirName = '1d_corrosion/1d_corrosion.csv' # directory containing csv file

params_data = readDatFile('dakota_tabular.dat')

params_data[0]+=  ["mass_loss","elapsed"]#['inital_total_cr',"final_total_cr","inital_metal_cr","final_metal_cr","corrosion_depth"]
#
for i in range(1,len(params_data)):
    cw_dir = workdir+str(i)+"/"
    data   = np.asarray(getRawData(cw_dir+csv_dirName,','))

    end_time = data[0,-1]
    if (end_time < 3.6e6):
        print("Incomplete simulation", (i - 1)/6)
    inital_metal_cr = data[1,1]
    final_metal_cr  = data[1,-1]
    mass_loss = final_metal_cr - inital_metal_cr
    params_data[i]+= [mass_loss,end_time]


x_count = len(params_data[0])-2
#
param_names = params_data[0]
sens_data   = np.asarray(params_data[1:])

X = sens_data[:,:x_count]
Y = sens_data[:,x_count:]

Sens_vals = []
for i in range(x_count):
    x = np.hstack([X[0,i],X[1+6*i:1+6*(i+1),i]])
    y = np.hstack([Y[0,0],Y[1+6*i:1+6*(i+1),0]])#mf-m0
    t = np.hstack([Y[0,1],Y[1+6*i:1+6*(i+1),1]])

    x = [x[i] for i,xi in enumerate(x) if t[i]==3600000 ]
    y = [y[i] for i,yi in enumerate(y) if t[i]==3600000 ]

    P = np.polyfit(x,y,1)
    S = P[0]*np.mean(x)/np.mean(y)
    Sens_vals+=[S]

Sens_vals = np.asarray(Sens_vals)
Sens_vals/=np.max(np.abs(Sens_vals))
print("Sensitivi = ",Sens_vals)


plt.close('all')
PATH = '/blue/michael.tonks/chaitanya.bhave/yellowjacket/yj_first_paper/Figures/'

#set plot dimensions
fig = plt.figure(figsize=(3.0,1.85),dpi=500)
plt.rcParams.update({'font.family':'Arial'})
plt.rc('font', family='serif',weight='bold')


ax = fig.add_subplot(111)

xlabels = ["$i_0$","$D_{Ni}^{metal}$","$D_{Ni}^{melt}$","$D_{Cr}^{metal}$","$D_{Cr}^{melt}$","$\delta_{int}$","$l$","$\sigma_{int}$"]
ax.bar(xlabels,Sens_vals,color='b' )

ymin=-0.1
ymax= 1.1
ax.set_ylim(ymin,ymax)

ax.set_xticks( xlabels )
ax.set_xticklabels( xlabels ,fontsize=6)

#
yticks = np.linspace(ymin,ymax,7)
ax.set_yticks( yticks)
ax.set_yticklabels( np.around(yticks,1),fontsize=6)

ax.tick_params(axis='x',direction='in',length=5)
ax.tick_params(axis='y',direction='in',which='major',length=5)

#Set axis spine widths
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)

ax.margins(x=0, y=0)
ax.set_ylabel('Normalized Sensitivities',fontsize=7,fontweight='bold')

legend_properties = {'weight':'bold','size':6}
lgd = ax.legend(["Mass loss"],  prop=legend_properties,framealpha=0)


# ax.set_xlabel('Parameter',fontsize=7,fontweight='bold')
fig.tight_layout(pad=0.3)  #Removes unnecessary padding from figure. Usually makes figure look much better
plt.savefig('1d_sensitivity.png',dpi=500,transparent=True)
plt.savefig(PATH+'1d_sensitivity.png',dpi=500,transparent=True)

"""

plt.ylim(-0.1,1.1)
plt.xticks(fontsize=12)  # arbitrary chosen
plt.yticks(np.linspace(-0.1,1.1, 7),fontsize=12)
ax = plt.axes()
ax.tick_params(axis='x',direction='in',pad=6)
ax.tick_params(axis='y',direction='in',pad=6)


plt.ylabel("Normalized Sensitivity",fontsize=15,fontweight='bold',labelpad=12)
plt.xlabel("Model parameters and material properties",fontsize=15,fontweight='bold',labelpad=12)

plt.savefig('1d_sensitivity.png',dpi=500,transparent=True)
 """

# plt.show()
