import matplotlib
import matplotlib.pyplot as plt
import csv
import math
import numpy as np
from scipy import stats

# close existing plots
plt.close('all')
PATH = '/blue/michael.tonks/chaitanya.bhave/yellowjacket/yj_first_paper/Figures/'
#set plot dimensions
fig = plt.figure(figsize=(3.0,1.85),dpi=500)
plt.rcParams.update({'font.family':'Arial'})
plt.rc('font', family='serif',weight='bold')

def getRawData(fileName, delim): # extract raw data from csv file
    rawData = []
    with open(fileName, 'r') as f:
        CSVReader = csv.reader(f, delimiter = delim, skipinitialspace = True)
        labels = next(CSVReader)
        array_size = len(labels)
        rawData = [ [] for i in range(array_size)]
        for row in CSVReader:
            for (i,val) in enumerate(row):
                try:
                    rawData[i].append(float(val) )
                except:
                    pass
    return (labels,rawData)
def get_plot_vars(training_data,idx):
    T = training_data[0]
    H = T/3600 #time in hours
    M = training_data[idx]
    dM = ((M[1] - M[1:]))
    return(H,dM)


## MASS LOSS PLOT
ax = fig.add_subplot(111)

labels,raw_data = getRawData('/blue/michael.tonks/chaitanya.bhave/yellowjacket/yj_first_paper/data_files/Ni5Cr_grain_size_raw.csv',',')
data = (raw_data)

D = np.asarray(data[0])
A_f = np.asarray(data[1])

# ax.plot(D,A_f,'r',linewidth=1.5)

mean_af = np.sum(A_f*D)
print("Mean grain size by area is = ",mean_af)

N_f = A_f/(D**2)
N_f = N_f/sum(N_f)
ax.plot(D,N_f,'b',linewidth=1.5)

mean_nf = np.sum(N_f*D)
err2 = (D - mean_nf)**2
sd = np.sqrt(np.sum(N_f*err2) )

print("Mean grain size by grain count is = ",mean_nf," SD = ",sd)


E = math.log(mean_nf**2/math.sqrt(mean_nf**2 + sd**2)) #2.095
SD = math.sqrt(math.log(1 + (sd/mean_nf)**2 ))#1.1697
D = np.linspace(0.1,210,100)
r = stats.norm(E,SD)
pdf = r.pdf(np.log(D) )
ax.plot(D, pdf,'b--',linewidth=1.5 )

mean_ln = math.exp(E)#math.exp(np.sum(pdf*np.log(D) )/np.sum(pdf))
print("Mean grain size by lognormal fit grain count is = ",mean_ln,np.sum(pdf) )

# ln_af = pdf*D**2
# print("Sum ln_af = ",sum(ln_af))
# ln_af = ln_af/sum(ln_af)
# ax.plot(D,ln_af,'r--',linewidth=1.5)

# legend_properties = {'weight':'bold','size':6}
# ax.legend(['A_f','N_f','lgN N_f'],prop=legend_properties)
# ax.set_xscale('log')
plt.show()
