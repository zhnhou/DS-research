import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from figure_spt_planck import *

def plot_field(fits_file, xra_deg, yra_deg, yticks=None, xticks=None, vrange=None, cbticks=None):

    mf = create_map_figure(fits_file)

    xra = [mf.ra0dec0[0]-xra_deg/2., mf.ra0dec0[0]+xra_deg/2.]
    yra = [mf.ra0dec0[1]-yra_deg/2., mf.ra0dec0[1]+yra_deg/2.]

    map2d = np.flipud(mf.cut_map(xra,yra)) * 1e6

    vmin = vrange[0]
    vmax = vrange[1]

    fig, ax = plt.subplots()

    im = plt.imshow(map2d, cmap=plt.get_cmap('bone'), extent=(xra[0],xra[1],yra[0],yra[1]), vmax=vmax, vmin=vmin,
                    interpolation='bicubic', aspect='equal')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.2)

    cbar = fig.colorbar(im, ticks=cbticks, orientation='vertical', cax=cax,
                        drawedges=False)

    cb_yticklabels = []
    for i in cbticks:
        cb_yticklabels.append("$"+str(i)+"$")

    cbar.solids.set_edgecolor("face")
    cbar.ax.set_ylabel(' ', fontsize=16)
    cbar.ax.axes.set_yticklabels(cb_yticklabels, fontsize=16)
    cbar.outline.set_edgecolor('black')
    cbar.outline.set_linewidth(1.2)

    #plt.setp(ax.xaxis.get_ticklines(), 'markersize', 3, 'markeredgewidth', 2)
    #plt.setp(ax.yaxis.get_ticklines(), 'markersize', 3, 'markeredgewidth', 2)

    ax.set_xticks(xticks)
    ax.set_yticks(yticks)

    if not (xticks is None):
        xticklabels = []
        for i in xticks:
            xticklabels.append("$"+str(i)+"^{\circ}$")
        ax.axes.set_xticklabels(xticklabels, fontsize=16)
        ax.set_xlabel(r"$\mathrm{RA}$", fontsize=16)

    if not (yticks is None):
        yticklabels = []
        for i in yticks:
            yticklabels.append("$"+str(i)+"^{\circ}$")
        ax.axes.set_yticklabels(yticklabels, fontsize=16)
        ax.set_ylabel(r'$\mathrm{Dec}$', fontsize=16)

    #fig_file = 'map_hfi143_filter.pdf'
    #plt.savefig(fig_file, format='pdf', transparent=True)
    #plt.clf()
