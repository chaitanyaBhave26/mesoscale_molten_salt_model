import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

fig, (ax1,ax2,ax3) = plt.subplots(3,figsize=(1.0, 2.885),dpi = 500,gridspec_kw={'height_ratios': [1, 1,1]})

cmap = mpl.cm.jet
norm = mpl.colors.Normalize(vmin=0, vmax=9)

cbar = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                norm=norm,
                                orientation='vertical')

cb_ticks = [0,3,6,9]
cbar.ax.tick_params(labelsize=6)
cbar.set_ticks(cb_ticks)
cbar.ax.set_yticklabels(cb_ticks )
cbar.ax.set_ylabel("Grain OPs",fontsize=8,fontweight='bold')

ax2.axis('off')#set_aspect('auto')

cmin=0.008
cmax = 0.054
cmap = mpl.cm.jet
norm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)

cbar = mpl.colorbar.ColorbarBase(ax3, cmap=cmap,
                                norm=norm,
                                orientation='vertical')

cb_ticks = np.round(np.linspace(cmin,cmax,4),3)#[0.008,3,6,0.054]
print(cb_ticks)
cbar.ax.tick_params(labelsize=6)
cbar.set_ticks(cb_ticks)
cbar.ax.set_yticklabels( ['0.8','2.3','3.9','5.4',])

cbar.ax.set_ylabel(r"$\bf{\times 10^{-2} \Delta c_{Cr}}$",fontsize=8,fontweight='bold')


#TIGHT LAYOUT ADJUSTS BORDERS AND PADDING TO GIVE BEST LOOKING IMAGE
plt.tight_layout()

#SAVE FIGURE WITH DPI=500 AND TRANSPARENT BACKGROUND
fig.savefig('/blue/michael.tonks/chaitanya.bhave/yellowjacket/yj_first_paper/Figures/3d_colorbar_combined.png',dpi=500,transparent=True)             #Remember to create the folder pyrender to store images in!!
plt.close()
