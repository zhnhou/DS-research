import numpy as np
import matplotlib.pyplot as plt
import healpy as hp

class Healpix_SPTField(object):
    def __init__(self, nside, background=None):
        self.nside = nside
        
        if not (background is None):
            self.background = background
        else:
            self.background = np.zeros(hp.nside2npix(self.nside))
            self.background[:] = -1.6375e30

    def vec2pix_ring(self, x, y, z):
        return hp.vec2pix(self.nside,x,y,z)

    def create_background(self):
#        self.projorth = hp.projector.OrthographicProj(xsize=1600, half_sky=True, rot=[0,270,0], coord=['G','C'])
#        self.projorth.set_proj_plane_info(xsize=1600, half_sky=True)
#        self.orth_map = self.projorth.projmap(self.background, self.vec2pix_ring)
#        plt.imsave('test.png', self.orth_map, cmap="RdBu_r", vmin=-400, vmax=400)

        hp.orthview(map=self.background, rot=[0,270,0], coord=['G','C'], half_sky=True, notext=True, min=-400, max=400, xsize=1600)
        hp.graticule()

        self.create_spt_area()       
  
    def create_spt_area(self):
        
        num_ra = 1000
        num_dec = 200

        ra = np.arange(num_ra, dtype=np.float32)/num_ra * (360.0+106.98 - 297.84) + 297.84
        ip = np.where(ra > 360.0)
        ra[ip] = ra[ip] - 360.00

        dec = 90.0 - (np.arange(num_dec, dtype=np.float32)/num_dec * (66.15 - 38.92) - 66.15)

        ra[:] = ra[:] * np.pi/180.0
        dec[:] = dec[:] * np.pi/180.0

        ra_list = []
        dec_list = []
        for i in np.arange(num_ra):
            ra_list.append( ra[i] )
            dec_list.append( dec[0] )

        for i in np.arange(num_dec):
            ra_list.append( ra[num_ra-1] )
            dec_list.append( dec[i] )

        for i in np.arange(num_ra):
            ra_list.append( ra[num_ra-1-i] )
            dec_list.append( dec[num_dec-1] )

        for i in np.arange(num_dec):
            ra_list.append( ra[0] )
            dec_list.append( dec[num_dec-1-i] )

        hp.projplot(dec_list, ra_list, 'k')
