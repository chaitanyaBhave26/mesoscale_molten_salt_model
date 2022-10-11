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


# user-defined parameters
workdir = 'workdir.'

csv_dirName = '1d_corrosion/1d_corrosion.csv' # directory containing csv file

params_data = readDatFile('dakota_tabular.dat')

params_data[0]+=  ["mass_loss","elapsed"]#['inital_total_cr',"final_total_cr","inital_metal_cr","final_metal_cr","corrosion_depth"]
#
for i in range(1,len(params_data)):
    cw_dir = workdir+str(i)+"/"
    data   = np.asarray(getRawData(cw_dir+csv_dirName,',')[1] )

    end_time = data[0,-1]
    if (end_time < 3.6e6):
        print("Incomplete simulation", (i - 1)/6)
    inital_metal_cr = data[1,1]
    final_metal_cr  = data[1,-1]
    mass_loss = final_metal_cr - inital_metal_cr
    params_data[i]+= [mass_loss,end_time]

plt.close('all')
PATH = '/blue/michael.tonks/chaitanya.bhave/yellowjacket/yj_first_paper/Figures/'
#set plot dimensions

plt.rcParams.update({'font.family':'Arial'})
plt.rc('font', family='sans-serif',weight='bold')

##Cr concentration vs mass loss
fig = plt.figure(figsize=(3.0,1.85),dpi=500)
ax = fig.add_subplot(111)


c_Cr_plot = [0.05,0.01,0.02,0.1,0.2]
D_Cr_plot = [3e-5,1e-5,2e-5,4e-5,5e-5]

param_names = params_data[0]
sens_data   = np.asarray(params_data[1:])

D_Cr = sens_data[:,0]
mu_Ni = sens_data[:,1]

Y = sens_data[:,-2]

Y_c = []
Y_D = []
for (i,y) in enumerate(Y):
    if mu_Ni[i] == mu_Ni[0]:
        Y_D += [-y]
    if D_Cr[i] == D_Cr[0]:
        Y_c += [-y]

#Trendline
sqrt_D = np.sqrt(D_Cr_plot)
D = np.linspace(0,1e-4,6)

density = 8.7e3 #mg/cm^3
V_a = 1.1131e-11*6.02214076e23*1e-12 #cm3/mol
density = 51.9961e3/V_a #mg/cm3
c0 = 0.05 - 0.00177

dM_analy = 2*c0*density*np.sqrt(D*1e-8*3.6e6/math.pi)
P2, = ax.plot( np.sqrt(D),dM_analy,'k-',linewidth=1.5,label="Analytical model") #,["",]

S1,C1 = np.polyfit(sqrt_D,Y_D,1)
S2,C2 = np.polyfit(np.sqrt(D),dM_analy,1)

err = (S2-S1)*100/S2
print("Error = ",err)
# P = np.polyfit( sqrt_D,Y_D,1)

y = Y_D
y_fit = 2*c0*density*np.sqrt( np.asarray(D_Cr_plot)*1e-8*3.6e6/math.pi)
SStot = sum((y-np.mean(y))**2);
SSres = sum((y-y_fit)**2);
Rsq = 1-SSres/SStot;
#
print(Rsq)
# P2, = ax.plot( np.sqrt(D),np.polyval(P,np.sqrt(D) ),'k--',linewidth=1.5,label="Trendline") #,["",]

#Actual plot
P1, = ax.plot(sqrt_D,Y_D,'bd',markersize=3,label="Time = 1000 hours")



# plt.figure()
# plt.plot(D_Cr_plot,Y_c,'r*')

xmin = 2e-3
xmax = 8e-3
ymin = 0
ymax = 0.8

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

xticks= np.linspace(xmin,xmax, 5)
ax.set_xticks( xticks )
ax.set_xticklabels( np.around(xticks*1e3,1) ,fontsize=6)

#
yticks = np.linspace(ymin,ymax, 5)
ax.set_yticks( yticks)
ax.set_yticklabels( np.around(yticks,1),fontsize=6)

ax.tick_params(axis='x',direction='in',length=5)
ax.tick_params(axis='y',direction='in',which='major',length=5)

#Set axis spine widths
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)

ax.set_ylabel('Mass loss (mg/cm$^2$)',fontsize=7,fontweight='bold')
ax.set_xlabel(r'$\sqrt{Cr \; diffusivity} ( \times 10^{-3} \sqrt{\mu m^2/s} ) $',fontsize=7,fontweight='bold')

legend_properties = {'weight':'bold','size':6}
lgd = ax.legend(handles=[P2,P1],prop=legend_properties,ncol=2,framealpha=0)

fig.tight_layout(pad=0.3)  #Removes unnecessary padding from figure. Usually makes figure look much better

fig.savefig(PATH+'1d_verification_fig4.png',dpi=500, transparent=True)
