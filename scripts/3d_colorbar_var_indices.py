import matplotlib.pyplot as plt
import matplotlib as mpl

fig, (ax,ax2) = plt.subplots(2,figsize=(1.0, 2),dpi = 500,gridspec_kw={'height_ratios': [2, 1]})

cmap = mpl.cm.jet
norm = mpl.colors.Normalize(vmin=0, vmax=9)

cbar = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                norm=norm,
                                orientation='vertical')

cb_ticks = [0,3,6,9]
cbar.ax.tick_params(labelsize=4)
cbar.set_ticks(cb_ticks)
cbar.ax.set_yticklabels(cb_ticks )
cbar.ax.set_ylabel("Grain OPs",fontsize=6,fontweight='bold')

ax2.axis('off')
x0 = 0
y0 = -1

#TIGHT LAYOUT ADJUSTS BORDERS AND PADDING TO GIVE BEST LOOKING IMAGE
plt.tight_layout()

#SAVE FIGURE WITH DPI=500 AND TRANSPARENT BACKGROUND
fig.savefig('/blue/michael.tonks/chaitanya.bhave/yellowjacket/yj_first_paper/Figures/3d_cb_var_indices.png',dpi=500)#,transparent=True)             #Remember to create the folder pyrender to store images in!!
plt.close()
