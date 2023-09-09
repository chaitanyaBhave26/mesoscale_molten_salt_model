from MultiExodusReader import MultiExodusReader

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PolyCollection
import matplotlib
import numpy as np
from time import time
import os
import math

t1 = time()
#EXODUS FILE FOR RENDERING
#ANY CHARACTER(S) CAN BE PLACED IN PLACE OF THE *, EG. 2D/grain_growth_2D_graintracker_out.e.1921.0000 or 2D/grain_growth_2D_graintracker_out.e-s001
filenames = 'workdir.'+str(100)+'/ni20_corr/ni20_corr.e*'

#READ EXODUS FILE SERIES WITH MultiExodusReader
MF = MultiExodusReader(filenames)

x,y,z,c = MF.get_data_at_time('c_Cr',MF.global_times[-1])               #Read coordinates and variable value --> Will be parallelized in future


X = x[:,0]
C = (c>np.max(c)*math.erf(1) )

sorted_list = (sorted(zip(X,C),key = lambda x: x[0] ) )
min_x =1000
for (x_val,c_val) in sorted_list:
    # x_val,c_val = point
    if x_val < min_x and not c_val:
        min_x = x_val

print(100,min_x)

"""
for i in range(1,101):
    filenames = 'workdir.'+str(i)+'/ni20_corr/ni20_corr.e*'

    #READ EXODUS FILE SERIES WITH MultiExodusReader
    MF = MultiExodusReader(filenames)

    x,y,z,c = MF.get_data_at_time('c_Cr',MF.global_times[-1])               #Read coordinates and variable value --> Will be parallelized in future


    X = x[:,0]
    C = (c>np.max(c)*math.erf(1) )

    sorted_list = (sorted(zip(X,C),key = lambda x: x[0] ) )
    min_x =1000
    for (x_val,c_val) in sorted_list:
        # x_val,c_val = point
        if x_val < min_x and not c_val:
            min_x = x_val

    print(i,min_x)

    # print(time()-t1)
 
    #GENERATE COORDINATES ARRAY THAT STORES X AND Y POINTS TOGETHER
    coords = np.asarray([ np.asarray([x_val,y_val]).T for (x_val,y_val) in zip(x,y) ])

    #GENERATE FIGURE WINDOW
    fig, ax = plt.subplots(figsize=(3,1) )

    #USE POLYCOLLECTION TO DRAW ALL POLYGONS DEFINED BY THE COORDINATES
    p = PolyCollection(coords, cmap=matplotlib.cm.coolwarm, alpha=1)      #Edge color can be set if you want to show mesh

    #COLOR THE POLYGONS WITH OUR VARIABLE
    p.set_array(np.array(c) )

    #ADD THE POLYGON COLLECTION TO AXIS --> THIS IS WHAT ACTUALLY PLOTS THE POLYGONS ON OUR WINDOW
    ax.add_collection(p)

    #Draw line showing corrosion depth
    ax.axvline(x=min_x, color='k', ymin=0, ymax=30, linestyle='--',linewidth=1.0)



    #FIGURE FORMATTING

    #SET X AND Y LIMITS FOR FIGURE --> CAN USE x,y ARRAYS BUT MANUALLY SETTING IS EASIER
    ax.set_xlim([0,250])
    ax.set_ylim([0,30])
    #SET ASPECT RATIO TO EQUAL TO ENSURE IMAGE HAS SAME ASPECT RATIO AS ACTUAL MESH
    ax.set_aspect('equal')

    #SET XTICKS AND YTICKS
    xticks = np.linspace(0,250, 6)
    ax.set_xticks( xticks)
    ax.set_xticklabels( xticks.astype(int),fontsize=6)


    yticks = np.linspace(0,30, 3)
    ax.set_yticks( yticks)
    ax.set_yticklabels( np.around(yticks,1),fontsize=6)

    #ADD A COLORBAR, VALUE SET USING OUR COLORED POLYGON COLLECTION
    fig.colorbar(p,label="$c_{Cr}$",ticks=[np.amin(c),np.amax(c)*math.erf(1),np.amax(c)])

    #STORE FIGURE IN 2D FOLDER, AND THE NAME ENDS WITH THE INDEX OF THE RENDERED FRAME. DPI = 500 AND TRANSPARENT BACKGROUND
    fig.savefig('2d_render/2d_render_'+str(i)+'.png',dpi=500,transparent=True )

    #CLOSE FIGURE AFTER YOU ARE DONE WITH IT. OTHERWISE ALL GENERATED FIGURES WILL BE HELD IN MEMORY TILL SCRIPT FINISHES RUNNING
    plt.close(fig)

    del fig
    
    del MF
    """