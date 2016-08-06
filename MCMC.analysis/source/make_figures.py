import numpy as np
from scipy.special import *

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
        self.ylabel = r'$10^2r_s(z_{\mathrm{drag}}) / D_v(0.57)$'

        xticklabels = []
        for x in self.xticks:
            xticklabels.append(r"$"+str(x)+"$")

        yticklabels = []
        for y in self.yticks:
            yticklabels.append(r"$"+str(y)+"$")

        self.xticklabels = xticklabels
        self.yticklabels = yticklabels

class Figure_CMBBAOH0(object):
    def __init__(self, cmb_type, chain_path, z_param='Omegamh2', NoColorBar=False, NoYticks=False):

        self.chain_path = chain_path
        self.cmb_type = cmb_type
        self.x_param = 'H0'
        self.y_param = 'rsDv057'
        self.z_param = 'z_param'

        self.NoColorBar = NoColorBar
        self.NoYticks = NoYticks

    
