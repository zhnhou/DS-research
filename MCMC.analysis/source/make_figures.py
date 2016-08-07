import numpy as np
from scipy.special import *
from mcmc_analysis import *
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
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

    def __init__(self, BAO_mean=7.315, BAO_err=0.118, H0_mean=73.8, H0_err=2.4):
        self.BAO_mean = BAO_mean
        self.BAO_err = BAO_err
        self.H0_mean = H0_mean
        self.H0_err = H0_err

        self.xmin = 60
        self.xmax = 90
        self.ymin = 6.5
        self.ymax = 9.5

        self.xticks = [60, 70, 80, 90]
        self.yticks = [7, 8, 9]

        self.gray_contour = ['lightgray', 'gray']

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
    def __init__(self, cmb_type, chain_path, z_param='Omegamh2', NoColorBar=False, NoYticks=False, 
        NoBAO=False, NoH0=False, pdf_file='figure.pdf'):

        self.chain_path = chain_path
        self.cmb_type = cmb_type
        self.x_param = 'H0'
        self.y_param = 'rsDv057'
        self.z_param = z_param
        self.pdf_file = pdf_file

        self.NoColorBar = NoColorBar
        self.NoYticks = NoYticks
        self.NoBAO = NoBAO
        self.NoH0 = NoH0

        self.fp = FigProp_CMBBAOH0()

    def process_chain(self):
        
        feature = ['H0', 'Dv_0.57', 'rs_zdrag', self.z_param]
        self.mcmc = mcmc_analysis(chain_path=self.chain_path, feature=feature)
        
        # 1.022 is the correction for BAO observables
        self.mcmc.chain['rsDv057'] = 100.0 * self.mcmc.chain['rs_zdrag'] / self.mcmc.chain['Dv_0.57'] * 1.0220

    def add_2DGaussian_prior(self):
        nbins = 401

        i = 0 
        for err_scale in [2.4477, 1.51]:
            self.ax.add_patch(Ellipse([self.fp.H0_mean, self.fp.BAO_mean], 2*err_scale*self.fp.H0_err, 2*err_scale*self.fp.BAO_err, fill=True, color=self.fp.gray_contour[i], zorder=i))
            i += 1

        plt.plot([self.fp.xmin, self.fp.xmax], [self.fp.BAO_mean, self.fp.BAO_mean], zorder=2, color='k', linewidth=0.8)
        plt.plot([self.fp.xmin, self.fp.xmax], [self.fp.BAO_mean - self.fp.BAO_err, self.fp.BAO_mean - self.fp.BAO_err], zorder=2, linestyle='--', color='k', linewidth=0.8)
        plt.plot([self.fp.xmin, self.fp.xmax], [self.fp.BAO_mean + self.fp.BAO_err, self.fp.BAO_mean + self.fp.BAO_err], zorder=2, linestyle='--', color='k', linewidth=0.8)

        plt.plot([self.fp.H0_mean, self.fp.H0_mean], [self.fp.ymin, self.fp.ymax],  zorder=2, color='k', linewidth=0.8)
        plt.plot([self.fp.H0_mean - self.fp.H0_err, self.fp.H0_mean - self.fp.H0_err], [self.fp.ymin, self.fp.ymax], zorder=2, linestyle='--', color='k', linewidth=0.8)
        plt.plot([self.fp.H0_mean + self.fp.H0_err, self.fp.H0_mean + self.fp.H0_err], [self.fp.ymin, self.fp.ymax], zorder=2, linestyle='--', color='k', linewidth=0.8)


    def create_figure(self):

        fig, ax = plt.subplots()
        
        self.ax = ax

        ax.set_position([0.1,0.125,0.8*13.0/16,0.8])

        x_select, y_select, z_select = self.mcmc.create_2d_scatter(self.x_param, self.y_param, self.z_param, num_select=10000)
        scatter = plt.scatter(x_select, y_select, c=z_select, s=1, vmin=0.105, vmax=0.155, edgecolors='face', cmap=plt.get_cmap('YlGnBu_r'), zorder=3)

        X, Y, Z = self.mcmc.create_2d_gmm(self.x_param, self.y_param, do_weight=False)

        z, clevels = calc_clevel_gmm(Z)
        plt.contour(X, Y, z, levels=clevels, colors='k', zorder=4)


        if (not (self.NoBAO or self.NoH0)):
            self.add_2DGaussian_prior()

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

        ax.set_xlim([self.fp.xmin, self.fp.xmax])
        ax.set_ylim([self.fp.ymin, self.fp.ymax])

        ax.set_xticks(self.fp.xticks)
        ax.set_yticks(self.fp.yticks)
    
        ax.set_xlabel(self.fp.xlabel, fontsize=16)
        ax.set_ylabel(self.fp.ylabel, fontsize=16)

        ax.set_xticklabels(self.fp.xticklabels, fontsize=16)
        ax.set_yticklabels(self.fp.yticklabels, fontsize=16) 

        plt.savefig(self.pdf_file, format='pdf')
