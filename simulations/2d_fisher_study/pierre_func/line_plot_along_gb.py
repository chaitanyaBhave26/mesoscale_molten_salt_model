#Plots concentration profile along GB
from MultiExodusReader import MultiExodusReader

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PolyCollection
import matplotlib
import numpy as np
from time import time
import os
import math
def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]

def line_plot_data(filenames,P1,P2,n_points):
    #CREATE MULTIEXODUSREADER OBJECT
    MF = MultiExodusReader(filenames)
    x,y,z,c = MF.get_data_at_time('c',MF.global_times[-1])           #Read coordinates and variable value --> Will be parallelized in future

    #GET CENTER VALUE FOR EACH CELL
    X_mean = np.mean(x,1)
    Y_mean = np.mean(y,1)

                                                  #NUM OF POINTS TO SAMPLE ALONG LINE

    #GET DISTANCE BETWEEN END POINTS
    D = math.sqrt(sum([ (P2[i] - P1[i])**2 for i in range(len(P1)) ]) )
    #GET DISTANCES FOR EACH SAMPLE POINT FROM P1
    d = np.linspace(0,D,n_points)

    #GET X AND Y CO-ORDINATES FOR SAMPLE POINTS
    P_x = P1[0] + d*(P2[0] -P1[0] )/D
    P_y = P1[1] + d*(P2[1] -P1[1] )/D

    #STORE X,Y COORDINATES IN AN ARRAY OF POINTS
    P = np.asarray([P_x,P_y]).T

    #FOR EACH SAMPLE POINT, FIND VARIABLE VALUE ON THE CLOSEST CELL
    plot_val = []
    for point in P:
        dst_sq = (X_mean-point[0])**2 + (Y_mean - point[1])**2
        nearest_idx = np.argmin(dst_sq)
        plot_val+=[c[nearest_idx] ]
    return (d,plot_val)


#FIGURE PARAMS
fig, ax = plt.subplots(1,1,figsize=(3,3),dpi = 500)

#PARAMETERS OF THE LINE WE WANT TO PLOT VALUES ALONG
P1 = [150-0,6]                                                          #END POINT
P2 = [0,6]                                                    #END POINT
n_points = (P1[0]-P2[0])

#workdir1
filenames = 'workdir.1/2_grain/2_grain.e'                #Star represents all files following this template
d,plot_val1 = line_plot_data(filenames,P1,P2,n_points)
#PLOT DISTANCE OF POINT FROM P1 (X-AXIS) VS VARIABLE VALUE (Y-AXIS)
ax.plot(d,plot_val1,'r')
er = np.abs(np.asarray(plot_val1)-math.erf(1))
idx = np.argmin(er)
print(d[idx])

#workdir2
filenames = 'workdir.2/2_grain/2_grain.e'                #Star represents all files following this template
d,plot_val2 = line_plot_data(filenames,P1,P2,n_points)
#PLOT DISTANCE OF POINT FROM P1 (X-AXIS) VS VARIABLE VALUE (Y-AXIS)
ax.plot(d,plot_val2,'g--')

#workdir3
filenames = 'workdir.3/2_grain/2_grain.e'                #Star represents all files following this template
d,plot_val3 = line_plot_data(filenames,P1,P2,n_points)
#PLOT DISTANCE OF POINT FROM P1 (X-AXIS) VS VARIABLE VALUE (Y-AXIS)
ax.plot(d,plot_val3,'b-.')

#workdir3
filenames = 'workdir.4/2_grain/2_grain.e'                #Star represents all files following this template
d,plot_val4 = line_plot_data(filenames,P1,P2,n_points)
#PLOT DISTANCE OF POINT FROM P1 (X-AXIS) VS VARIABLE VALUE (Y-AXIS)
ax.plot(d,plot_val4,'m-')

#workdir3
filenames = 'workdir.5/2_grain/2_grain.e'                #Star represents all files following this template
d,plot_val5 = line_plot_data(filenames,P1,P2,n_points)
#PLOT DISTANCE OF POINT FROM P1 (X-AXIS) VS VARIABLE VALUE (Y-AXIS)
ax.plot(d,plot_val5,'c-.')

# err1 = (plot_val5-plot_val1)/

# #Fisher prediction
# y = d
# t = 1e6
# fisher_plot_c_y=np.exp(-(1.33*y)/((2.07e-7*t)**(1/2)*(1750.0/(2.07e-7*t)**(1/2))**(1/2)))
# ax.plot(d,fisher_plot_c_y,'k--')

#PLOT FORMATTING
xmin = 0
xmax = P1[0]
ymin = 0.0
ymax = 1.0


ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])

xticks=[ int(tick) for tick in np.linspace(0,xmax, 6)]
ax.set_xticks( xticks )
ax.set_xticklabels( xticks ,fontsize=6,fontweight='bold')
yticks = np.around(np.linspace(ymin,ymax, 7),1)
ax.set_yticks( yticks)
ax.set_yticklabels( yticks,fontsize=6,fontweight='bold')
ax.tick_params(axis='x',direction='in',length=5)
ax.tick_params(axis='y',direction='in',which='major',length=5)
ax.margins(x=0, y=0)
ax.set_ylabel('c',fontsize=7,fontweight='bold')
ax.set_xlabel(r'Distance ($\mu$m)',fontsize=7,fontweight='bold')


#TIGHT LAYOUT ADJUSTS BORDERS AND PADDING TO GIVE BEST LOOKING IMAGE
plt.tight_layout()

#SAVE FIGURE WITH DPI=500 AND TRANSPARENT BACKGROUND
fig.savefig('2d_lineplot.png',dpi=500,transparent=True)             #Remember to create the folder pyrender to store images in!!
plt.close()
