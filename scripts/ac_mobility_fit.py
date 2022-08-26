import matplotlib.pyplot as plt
import numpy as np
from matplotlib.legend_handler import HandlerTuple

# close existing plots
plt.close('all')
PATH = '/blue/michael.tonks/chaitanya.bhave/yellowjacket/yj_first_paper/Figures/'
#set plot dimensions
fig = plt.figure(figsize=(3.0,1.85),dpi=500)
ax = fig.add_subplot(1,1,1)
plt.rcParams.update({'font.family':'Arial'})
plt.rc('font', family='sans-serif',weight='bold')

dx = 0.5
X = np.linspace(-3,3,500)

int_width = 0.5
eta = 0.5*(1 - np.tanh( 2*X/int_width ) )
P1,=ax.plot(X,eta,'r-',linewidth=1.5)

eta = 0.5*(1 - np.tanh( 2*(X+dx)/int_width ) )
P2,=ax.plot(X,eta,'r--',linewidth=1.5)

xmin = -1.0
xmax = 1.0
ymin = -0.1
ymax = 1.1
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

xticks=[ int(tick) for tick in np.linspace(xmin,xmax, 7)]
ax.set_xticks(xticks)  # arbitrary chosen
ax.set_xticklabels(xticks,fontsize=6,fontweight='bold')

yticks = np.around(np.linspace(ymin,ymax,7),1)
ax.set_yticks(yticks)
ax.set_yticklabels(yticks,fontsize=6,fontweight='bold')

ax.tick_params(axis='x',direction='in',length=5)#,pad=6,)
ax.tick_params(axis='y',direction='in',length=5)#,pad=6)

legend_properties = {'weight':'bold','size':6}
legs = ["Sharp interface (t=0)","Sharp interface (t=dt)","Phase-field (t=0)","Phase-field (t=dt)"]
lgd = plt.legend( ['t=0','t = dt'],frameon=False, prop=legend_properties, handlelength=3, handler_map={tuple: HandlerTuple(ndivide=None)})

plt.margins(x=0, y=0)
plt.xlabel("X-axis ($\mu$m)",fontsize=6,fontweight='bold')#,labelpad=15)
plt.ylabel("Metal phase",fontsize=6,fontweight='bold')#,labelpad=15)

#Set axis spine widths
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1.0)

fig.tight_layout()
fig.savefig(PATH+'AC_mobility_fit.png',dpi=500, transparent=True)
