import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

fig, (ax,ax2) = plt.subplots(2,figsize=(1.0, 2),dpi = 500,gridspec_kw={'height_ratios': [2, 1]})
cmin=0.008
cmax = 0.054
cmap = mpl.cm.jet
norm = mpl.colors.Normalize(vmin=cmin, vmax=cmax)

cbar = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                norm=norm,
                                orientation='vertical')

cb_ticks = np.round(np.linspace(cmin,cmax,4),3)#[0.008,3,6,0.054]
print(cb_ticks)
cbar.ax.tick_params(labelsize=4)
cbar.set_ticks(cb_ticks)
cbar.ax.set_yticklabels( ['0.8','2.3','3.9','5.4',])

cbar.ax.set_ylabel(r"$\bf{\times 10^{-2} \Delta c_{Cr}}$",fontsize=6,fontweight='bold')

ax2.axis('off')#set_aspect('auto')

#TIGHT LAYOUT ADJUSTS BORDERS AND PADDING TO GIVE BEST LOOKING IMAGE
plt.tight_layout()

#SAVE FIGURE WITH DPI=500 AND TRANSPARENT BACKGROUND
fig.savefig('/blue/michael.tonks/chaitanya.bhave/yellowjacket/yj_first_paper/Figures/3d_cb_del_Cr.png',dpi=500)#,transparent=True)             #Remember to create the folder pyrender to store images in!!
plt.close()
