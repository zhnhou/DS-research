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


    def create_background(self):
        hp.orthview(map=self.background, rot=[0,270,0], coord=['G','C'])
