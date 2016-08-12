import numpy as np
from scipy.special import *
from mcmc_analysis import *
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.gridspec as gridspec
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


class Figure_ChainAnalysis(object):
    def __init__(self, chain_path, chain_prefix=None, feature=None, is_cosmomc=False, num_chain=8):
        self.mcmc = mcmc_analysis(chain_path=chain_path, chain_prefix=chain_prefix, feature=feature, is_cosmomc=is_cosmomc, num_chain=num_chain)

    def create_chain_burnin(self, parameter=None, ncol=4, nrow=2, wspace=0.0, hspace=0.0, left=0.1, right=0.95, num_step=2000, rescale_parameter=False):

        self.rescale_parameter = rescale_parameter

        if not (parameter is None):
            self.parameter = parameter
        else:
            self.parameter = mcmc.parameter
        gs = gridspec.GridSpec(nrow,ncol)
        gs.update(left=left, right=right, wspace=wspace, hspace=hspace)
        
        i = 0
        for irow in np.arange(nrow):
            for icol in np.arange(ncol):
                param = self.parameter[i]

                ax = plt.subplot(gs[irow,icol])
                

                x = self.mcmc.chain_original['step'][0:num_step]
                if self.rescale_parameter:
                    y = self.mcmc.chain_original[param][0:num_step]
                else:
                    mean = np.mean(self.mcmc.chain[param])
                    stddev = np.std(self.mcmc.chain[param])

                    y = (self.mcmc.chain_original[param][0:num_step] - mean) / stddev

#                pmax = np.amax(y)
#                pmin = np.amin(y)

                ax.plot(x, y)

                ax.set_xlim([0,num_step])
                ax.set_ylim([-4.5,4.5])

                if irow != nrow-1:
                    ax.set_xticklabels([])
                if icol != 0:
                    ax.set_yticklabels([])
                i += 1

        plt.show()




class FigProp_CMBBAOH0(object):

    def __init__(self, BAO_mean=7.315, BAO_err=0.118, H0_mean=73.8, H0_err=2.4):
        self.BAO_mean = BAO_mean
        self.BAO_err = BAO_err
        self.H0_mean = H0_mean
        self.H0_err = H0_err

        self.xmin = 61
        self.xmax = 89
        self.ymin = 6.5
        self.ymax = 9.5
        self.zmin = 0.105
        self.zmax = 0.155

        self.xticks = [65, 70, 75, 80, 85]
        self.yticks = [7, 8, 9]
        self.cbticks = [0.11,0.12,0.13,0.14,0.15]

        self.gray_contour = ['lightgray', 'gray']

        self.xlabel = r'$H_0\;[\mathrm{km\,s^{-1}\,Mpc^{-1}}]$'
        self.ylabel = r'$10^2r_s(z_{\mathrm{drag}}) / D_{\mathrm{v}}(0.57)$'
        self.cblabel = r'$\Omega_m h^2$'

        xticklabels = []
        for x in self.xticks:
            xticklabels.append(r"$"+str(x)+"$")

        yticklabels = []
        for y in self.yticks:
            yticklabels.append(r"$"+str(y)+"$")

        cbticklabels = []
        for z in self.cbticks:
            cbticklabels.append(r"$"+str(z)+"$")

        self.xticklabels = xticklabels
        self.yticklabels = yticklabels
        self.cbticklabels = cbticklabels

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
            self.ax.add_patch(Ellipse([self.fp.H0_mean, self.fp.BAO_mean], 2*err_scale*self.fp.H0_err, 2*err_scale*self.fp.BAO_err, fill=True, color=self.fp.gray_contour[i], zorder=i, linewidth=0))
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

        ax.set_position([0.1,0.125,0.850*13.0/16,0.850])

        x_select, y_select, z_select = self.mcmc.create_2d_scatter(self.x_param, self.y_param, self.z_param, num_select=10000)
        scatter = plt.scatter(x_select, y_select, c=z_select, s=1, vmin=self.fp.zmin, vmax=self.fp.zmax, edgecolors='face', cmap=plt.get_cmap('YlGnBu_r'), zorder=3)

        X, Y, Z = self.mcmc.create_2d_gmm(self.x_param, self.y_param, do_weight=False)

        z, clevels = calc_clevel_gmm(Z)
        plt.contour(X, Y, z, levels=clevels, colors='k', zorder=4)


        if (not (self.NoBAO or self.NoH0)):
            self.add_2DGaussian_prior()

        if (not self.NoColorBar):
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.2)
            
#            cbar = fig.colorbar(scatter, ticks=(-100,-50,0,50,100), orientation='vertical', cax=cax,
            cbar = fig.colorbar(scatter, ticks=self.fp.cbticks, orientation='vertical', cax=cax,
                   drawedges=False)

            cbar.solids.set_edgecolor("face")
            cbar.ax.set_ylabel(self.fp.cblabel, fontsize=18)
            cbar.ax.axes.set_yticklabels(self.fp.cbticklabels, fontsize=16)
            cbar.outline.set_edgecolor('black')
            cbar.outline.set_linewidth(1.2)

        ax.set_xlim([self.fp.xmin, self.fp.xmax])
        ax.set_ylim([self.fp.ymin, self.fp.ymax])

        ax.set_xticks(self.fp.xticks)
        ax.set_yticks(self.fp.yticks)
    
        ax.set_xlabel(self.fp.xlabel, fontsize=18)
        ax.set_ylabel(self.fp.ylabel, fontsize=18)

        ax.set_xticklabels(self.fp.xticklabels, fontsize=16)
        ax.set_yticklabels(self.fp.yticklabels, fontsize=16) 

        plt.savefig(self.pdf_file, format='pdf')
        os.system("convert -density 400 "+self.pdf_file+" "+self.pdf_file[:-4]+".png")


class Figure_CMBFeature(object):
    def __init__(self):
        self.wmap_
