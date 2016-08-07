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
#        self.projorth = hp.projector.SphericalProj()
#        self.orth_map = self.projorth.projmap(self.background, self.vec2pix_ring)
#        plt.imsave('test.png', self.orth_map, cmap="RdBu", vmin=-400, vmax=400)

        hp.orthview(map=self.background, rot=[0,270,0], coord=['G','C'], half_sky=True, notext=True, min=-400, max=400, xsize=1600)
        hp.graticule()
