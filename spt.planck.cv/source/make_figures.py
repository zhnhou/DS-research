import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import healpy as hp

def planck_color_map(self):

    cmap = ListedColormap(np.loadtxt("Planck_FreqMap_RGB.txt")/255.)
    cmap.set_bad("darkgray")
    cmap.set_under("white")

    return cmap

class Healpix_SPTField(object):
    def __init__(self, nside, background=None, input_coord='G', bad_pixels=None, offset=0, factor=1.0, sub=None, cbar=True, title=' '):
        self.nside = nside
        self.input_coord = input_coord
        self.bad_pixels = bad_pixels
        self.offset = offset
        self.factor = factor
        self.sub = sub
        self.cbar = cbar
        self.title = title
        
        if not (background is None):
            self.background = background
        else:
            self.background = np.zeros(hp.nside2npix(self.nside))
            self.background[:] = -1.6375e30

    def vec2pix_ring(self, x, y, z):
        return hp.vec2pix(self.nside,x,y,z)

    def create_background(self):

#        hp.orthview(map=self.background, rot=[0.0,270.0,0.0], coord=['G','C'], half_sky=True, notext=True, min=-400, max=400, xsize=3000)

        m = self.background * self.factor + self.offset * self.factor
        msinh = np.arcsinh(m * 0.5)/np.log(10.0)

        cbticks = np.array([-500.0,0.0,100.0,1000.0, 1e6])
        cbticklabels = np.array([r'$-500$',r'$0$',r'$100$',r'$1000$',r'$10^6$'])
        cbticks = np.arcsinh( (cbticks-80.0) * 0.5) / np.log(10.0)
         
        if not (self.bad_pixels is None):
            msinh[self.bad_pixels] = -1.6375e30
        
        if (self.cbar):
            orth = hp.orthview(map=msinh, rot=[0.0,270.0,0.0], coord=[self.input_coord,'C'], half_sky=True, notext=True, xsize=3000, cmap=self.planck_color_map(), min=-3.1, max=7, cbticks=cbticks, cbticklabels=cbticklabels, title=self.title, cbtitle=r'$\mathrm{\mu K}$', sub=self.sub, return_projected_map=True)
        else:
            orth = hp.orthview(map=msinh, rot=[0.0,270.0,0.0], coord=[self.input_coord,'C'], half_sky=True, notext=True, xsize=3000, cmap=self.planck_color_map(), min=-3.1, max=7, cbar=False, title=self.title, sub=self.sub, return_projected_map=True)

#        self.create_spt_area()
#        self.create_graticule([160,140,120,100], np.arange(0,360,45))
#        plt.savefig('test.png', format='png', dpi=600)

        return orth

    def save_figure(fig_file, dpi=600):
        form = fig_file[-3:]
        plt.savefig(fig_file, format=form, dpi=dpi)

    def create_graticule(self, lat, lon):

        for la in lat:
            theta = np.zeros(1000)
            theta[:] = la * np.pi/180
            phi = np.linspace(0,np.pi*2,num=1000)
            hp.projplot(theta, phi, 'gray', linewidth=0.3)
            hp.projplot(theta, phi, 'k', linewidth=1, linestyle='dotted')
        for lo in lon:
            theta = np.linspace(0,np.pi,num=1000)
            phi = np.zeros(1000)
            phi[:] = lo * np.pi/180
            hp.projplot(theta, phi, 'gray', linewidth=0.3)
            hp.projplot(theta, phi, 'k', linewidth=1, linestyle='dotted')
            
  
    def create_spt_area(self):
        
        num_ra = 1000
        num_dec = 200

        ra = np.arange(num_ra, dtype=np.float32)/num_ra * (360.0+105 - 300) + 300.0
        ip = np.where(ra > 360.0)
        ra[ip] = ra[ip] - 360.00

        dec = 90.0 - (np.arange(num_dec, dtype=np.float32)/num_dec * (65.0 - 40.00) - 65.0)

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


    def planck_color_map(self):

        '''

        Itab = [   0,  13,  26,  39,  52,  65,  76,  77,  88, 
                   101, 114, 127,   140, 153, 166, 179,   192, 
                   205, 218, 231, 255]
        Rtab = [   0,  10,  30,  80, 191, 228, 241, 241, 245, 
                   248, 249.9, 242.25,204, 165, 114, 127.5, 178.5, 
                   204, 229.5, 242.25, 252.45]
        Gtab = [   0,  20, 184, 235, 239, 240, 241, 241, 240, 
                   235, 204,   153,  76.5,  32,   0, 127.5, 178.5, 
                   204, 229.5, 242.25, 252.45]
        Btab = [ 255, 255, 255, 255, 250, 245, 212, 212, 175, 
                 130, 38.25, 12.75,  0,   32,  32, 153,   204,   
                 229.5, 242.25, 249.9, 255]

        ncolors = 256
        ii = np.arange(ncolors, dtype=np.float32)

        R = np.interp(ii, Itab, Rtab)
        G = np.interp(ii, Itab, Rtab)
        B = np.interp(ii, Itab, Rtab)

        '''

        cmap = ListedColormap(np.loadtxt("Planck_FreqMap_RGB.txt")/255.)
        cmap.set_bad("darkgray")
        cmap.set_under("white")

        return cmap


class Healpix_MollSPT(object):

    def __init__(self, background, vmin=None, vmax=None, cbar=True, sub=None, title=None, margins=None):

        self.background = background
        self.nside = hp.npix2nside(np.shape(background)[0])
        self.vmin = vmin
        self.vmax = vmax
        self.cbar = cbar
        self.sub = sub
        self.margins = margins

        if (title is None):
            self.title = ' '
        else:
            self.title = title

    def create_background(self):
        view = hp.mollview(map=self.background, xsize=3000, min=self.vmin, max=self.vmax, cbar=self.cbar, sub=self.sub, title=self.title, return_projected_map=True, margins=self.margins)

        return view


    def create_spt_area(self):
        
        num_ra = 1000
        num_dec = 200

        ra = np.arange(num_ra, dtype=np.float32)/num_ra * (360.0+105 - 300) + 300.0
        ip = np.where(ra > 360.0)
        ra[ip] = ra[ip] - 360.00

        dec = 90.0 - (np.arange(num_dec, dtype=np.float32)/num_dec * (65.0 - 40.00) - 65.0)

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

        rot = hp.Rotator(coord=['C','G'])
        theta, phi = rot(dec_list, ra_list)

        lines1 = hp.projplot(theta, phi, 'k')
        
        return lines1

