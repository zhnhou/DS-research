import numpy as np
from sptsz_map import *

class create_map_figure(object):
    def __init__(self, fits_file):
        self.fits_file = fits_file
        m = read_sptsz_fits(fits_file)

        self.ra0dec0 = [m['ra0'], m['dec0']]
        self.nside   = [m['nsidex'], m['nsidey']]
        self.map_data = np.flipud(m['map_data'])
        self.reso_arcmin = m['reso_arcmin']

        self.setup_coord()
        
    def setup_coord(self):
        reso_deg = self.reso_arcmin / 60.00

        self.xra = (self.ra0dec0[0] - 0.5*self.nside[0]*reso_deg, self.ra0dec0[0] + 0.5*self.nside[0]*reso_deg)
        self.yra = (self.ra0dec0[1] - 0.5*self.nside[1]*reso_deg, self.ra0dec0[1] + 0.5*self.nside[1]*reso_deg)

    def cut_map(self, xra, yra, replace=True):
        reso_deg = self.reso_arcmin / 60.00

        xarr = np.arange(0,self.nside[0])*reso_deg + self.ra0dec0[0] - 0.5*self.nside[0]*reso_deg
        yarr = np.arange(0,self.nside[1])*reso_deg + self.ra0dec0[1] - 0.5*self.nside[1]*reso_deg

        ipx = np.where((xarr >= min(xra)) & (xarr <= max(xra)))
        ipy = np.where((yarr >= min(yra)) & (yarr <= max(yra)))

        map2d_cut = self.map_data[ipy[0].min():(ipy[0].max()+1),ipx[0].min():(ipx[0].max()+1)]

        if (replace):
            self.map_data = map2d_cut
            self.xra = xra
            self.yra = yra

        return map2d_cut


