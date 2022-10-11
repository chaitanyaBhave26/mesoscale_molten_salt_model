# import packages
import matplotlib.pyplot as plt
import matplotlib
import csv
import math
import numpy as np
from scipy import stats
from IPython.display import set_matplotlib_formats
# set_matplotlib_formats('retina', quality=100)
#
# font = {'family' : 'serif',
#         'size'   : 24}
#
# matplotlib.rc('font', **font)

def getRawData(fileName, delim): # extract raw data from csv file
    rawData = []
    with open(fileName, 'r') as f:
        CSVReader = csv.reader(f, delimiter = delim, skipinitialspace = True)
        labels = next(CSVReader)
        array_size = len(labels)
        rawData = [ [] for i in range(array_size)]
        for row in CSVReader:
            for (i,val) in enumerate(row):
                rawData[i].append(float(val) )
    return (labels,rawData)


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
                # print("err")
                pass
        rawData+=[line]
    return rawData

plt.close('all')
PATH = '/blue/michael.tonks/chaitanya.bhave/yellowjacket/yj_first_paper/Figures/'
#set plot dimensions

plt.rcParams.update({'font.family':'Arial'})
plt.rc('font', family='sans-serif',weight='bold')

##Time vs dM --> 1D model vs analytical expression
fig = plt.figure(figsize=(3.0,1.85),dpi=500)
ax = fig.add_subplot(111)




labels,raw_data = getRawData('workdir.1/1d_corrosion/1d_corrosion.csv',',')
data = np.asarray(raw_data)
H,dM=get_plot_vars(data,1)

##ANALYTICAL CURVE
c0 = 0.05-0.00177
D  = 3e-5
density = 8.7e3 #mg/cm^3
V_a = 1.1131e-11*6.02214076e23*1e-12 #cm3/mol
density = 51.9961e3/V_a #mg/cm3
dM_analy = 2*c0*density*np.sqrt(D*1e-8*H*3600/math.pi)
ax.plot(H, dM_analy,'k',linewidth=1.5 )

##YELLOWJACKET RESULT
n = 40
ax.plot(H[1::n],dM[::n],'rd',markersize=3)
#Force last point plot
ax.plot(H[-1],dM[-1],'bd', zorder=10, clip_on=False,markersize=3)

error = np.sqrt( np.mean( (dM_analy[1:]-dM)**2 ) )
print("RMS error = ",error,"mg/cm^2")

xmin = H[0]
xmax = H[-1]
ymin = dM[0]
ymax = 0.5

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

xticks=[ int(tick) for tick in np.linspace(xmin,xmax, 6)]
ax.set_xticks( xticks )
ax.set_xticklabels( xticks ,fontsize=6)

#
yticks = np.linspace(ymin,ymax, 6)
ax.set_yticks( yticks)
ax.set_yticklabels( np.around(yticks,1),fontsize=6)

ax.tick_params(axis='x',direction='in',length=5)
ax.tick_params(axis='y',direction='in',which='major',length=5)

#Set axis spine widths
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)


ax.margins(x=0, y=0)
ax.set_ylabel('Mass loss (mg/cm$^2$)',fontsize=7,fontweight='bold')
ax.set_xlabel('Time (hours)',fontsize=7,fontweight='bold')

legend_properties = {'weight':'bold','size':6}
lgd = ax.legend(["Analytical model","1D PF model"],prop=legend_properties,ncol=2,framealpha=0)


ax.margins(x=0, y=0)
fig.tight_layout(pad=0.3)  #Removes unnecessary padding from figure. Usually makes figure look much better

fig.savefig(PATH+'1d_verification_fig1.png',dpi=500, transparent=True)
