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

int_widths = [4.0,2.0,1.0,0.5,0.2] #,0.1
dM_Cr_list = []

for l_pf in int_widths:
    file_name = '1d_corrosion/1d_corrosion'+str(l_pf)+'.csv'
    #Cr idx is 4
    data = np.asarray(getRawData(file_name,','))
    dM_Cr = data[4,-1]
    T = data[0,-1]
    dM_Cr_list+=[-dM_Cr]

ax.plot(int_widths,dM_Cr_list,'k*')

P = np.polyfit(int_widths,dM_Cr_list,1)
sharp_val = np.polyval(P,0)
print("Mass loss at 0 micron interface = ",sharp_val)
print("Error for 2 micron interface = ", (sharp_val-dM_Cr_list[1])*100/sharp_val," %" )
print("Error for 1 micron interface = ", (sharp_val-dM_Cr_list[2])*100/sharp_val," %" )
print("Error for 0.5 micron interface = ", (sharp_val-dM_Cr_list[3])*100/sharp_val," %" )



xmin = 0
xmax = 5
ymin = 0.59
ymax = 0.62

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

xticks=np.linspace(xmin,xmax, 6) #[ int(tick) for tick in np.linspace(xmin,xmax, 4)]
ax.set_xticks( xticks )
ax.set_xticklabels( xticks ,fontsize=6)

# #
yticks = np.linspace(ymin,ymax, 4)
ax.set_yticks( yticks)
ax.set_yticklabels( yticks,fontsize=6)

ax.tick_params(axis='x',direction='in',length=5)
ax.tick_params(axis='y',direction='in',which='major',length=5)

# #Set axis spine widths
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)


ax.margins(x=0, y=0)
ax.set_ylabel(r'Mass loss (mg/cm$^2$)',fontsize=7,fontweight='bold')
ax.set_xlabel(r'Interface width $l\ (\mu$m)',fontsize=7,fontweight='bold')

fig.tight_layout(pad=0.3)  #Removes unnecessary padding from figure. Usually makes figure look much better

fig.savefig(PATH+'1d_interface_width_mass_loss.png',dpi=500, transparent=True)

plt.show()
