import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

fig, (ax,ax2) = plt.subplots(2,figsize=(1.0, 2),dpi = 500,gridspec_kw={'height_ratios': [2, 1]})
# fig.subplots_adjust(bottom=0.5)
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

# cbar.ax.set_yticklabels( [r'$8.0\times10^{-3}$',r'$2.3\times10^{-2}$',r'$3.9\times10^{-2}$',r'$5.4\times10^{-3}$',] ,fontweight='bold')
cbar.ax.set_ylabel(r"$\bf{\times 10^{-2} \Delta c_{Cr}}$",fontsize=6,fontweight='bold')

# ax.set_aspect('auto')
ax2.axis('off')#set_aspect('auto')
# ax.annotate('', xy=(0, -0.1), xycoords='axes fraction', xytext=(1, -0.1),
# arrowprops=dict(arrowstyle="->", color='b'))
#
# x0 = 0
# y0 = -1
# ax2.annotate('', xy=(x0,y0 ), xycoords='axes fraction', xytext=(x0, y0+0.5),
# arrowprops=dict(arrowstyle="<-", color='b'))

#TIGHT LAYOUT ADJUSTS BORDERS AND PADDING TO GIVE BEST LOOKING IMAGE
plt.tight_layout()

#SAVE FIGURE WITH DPI=500 AND TRANSPARENT BACKGROUND
fig.savefig('/blue/michael.tonks/chaitanya.bhave/yellowjacket/yj_first_paper/Figures/3d_cb_del_Cr.png',dpi=500)#,transparent=True)             #Remember to create the folder pyrender to store images in!!
plt.close()
