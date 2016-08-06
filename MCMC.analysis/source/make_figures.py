import numpy as np
from scipy.special import *
from mcmc_analysis import *
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

def calc_clevel_gmm(logprob):
    z = np.exp(-logprob)
    nm = np.sum(z)
    z /= nm

    psort = np.sort(z.flatten())[::-1]
    pcum  = np.cumsum(psort)

    clevels = []
    sigma = [2.0, 1.0]
    ptes = erf(sigma / np.sqrt(2.0))

    for p in ptes:
        ind = np.where(pcum > p)[0][0]
        clevels.append(psort[ind])

    return z, clevels


class FigProp_CMBBAOH0(object):

    def __init__(self):
        self.xmin = 60
        self.xmax = 90
        self.ymin = 6.5
        self.ymax = 9.5

        self.xticks = [60, 70, 80, 90]
        self.yticks = [7, 8, 9]

        self.xlabel = r'$H_0\;[\mathrm{km\,s^{-1}\,Mpc^{-1}}]$'
        self.ylabel = r'$10^2r_s(z_{\mathrm{drag}}) / D_{\mathrm{v}}(0.57)$'

        xticklabels = []
        for x in self.xticks:
            xticklabels.append(r"$"+str(x)+"$")

        yticklabels = []
        for y in self.yticks:
            yticklabels.append(r"$"+str(y)+"$")

        self.xticklabels = xticklabels
        self.yticklabels = yticklabels

class Figure_CMBBAOH0(object):
    def __init__(self, cmb_type, chain_path, z_param='Omegamh2', NoColorBar=False, NoYticks=False, pdf_file='figure.pdf'):

        self.chain_path = chain_path
        self.cmb_type = cmb_type
        self.x_param = 'H0'
        self.y_param = 'rsDv057'
        self.z_param = z_param
        self.pdf_file = pdf_file

        self.NoColorBar = NoColorBar
        self.NoYticks = NoYticks

    def process_chain(self):
        
        feature = ['H0', 'Dv_0.57', 'rs_zdrag', self.z_param]
        self.mcmc = mcmc_analysis(chain_path=self.chain_path, feature=feature)
        self.mcmc.chain['rsDv057'] = 100.0 * self.mcmc.chain['rs_zdrag'] / self.mcmc.chain['Dv_0.57']


    def create_figure(self):

        fp = FigProp_CMBBAOH0()

        fig, ax = plt.subplots()

        ax.set_position([0.1,0.125,0.8,0.8])

        x_select, y_select, z_select = self.mcmc.create_2d_scatter(self.x_param, self.y_param, self.z_param, num_select=10000)
        scatter = plt.scatter(x_select, y_select, c=z_select, s=1, vmin=0.105, vmax=0.155, edgecolors='face', cmap=plt.get_cmap('YlGnBu_r'))

        X, Y, Z = self.mcmc.create_2d_gmm(self.x_param, self.y_param)

        z, clevels = calc_clevel_gmm(Z)

        plt.contour(X, Y, z, levels=clevels, colors='k')

#        ax.set_aspect((fp.xmax-fp.xmin)/(fp.ymax-fp.ymin))

        if (not self.NoColorBar):
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.2)
            
#            cbar = fig.colorbar(scatter, ticks=(-100,-50,0,50,100), orientation='vertical', cax=cax,
            cbar = fig.colorbar(scatter, orientation='vertical', cax=cax,
                   drawedges=False)

            cbar = fig.colorbar(scatter, orientation='vertical', cax=cax, drawedges=False)
            
            cbar.solids.set_edgecolor("face")
#            cbar.ax.set_ylabel(' ', fontsize=16)
#            cbar.ax.axes.set_yticklabels(["$-100$","$-50$","$0$","$50$","$100$"], fontsize=16)
            cbar.outline.set_edgecolor('black')
            cbar.outline.set_linewidth(1.2)

        ax.set_xlim([fp.xmin, fp.xmax])
        ax.set_ylim([fp.ymin, fp.ymax])

        ax.set_xticks(fp.xticks)
        ax.set_yticks(fp.yticks)
    
        ax.set_xlabel(fp.xlabel, fontsize=16)
        ax.set_ylabel(fp.ylabel, fontsize=16)

        ax.set_xticklabels(fp.xticklabels, fontsize=16)
        ax.set_yticklabels(fp.yticklabels, fontsize=16) 

        plt.savefig(self.pdf_file, format='pdf')
