import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import healpy as hp
from scipy.io.idl import readsav

def planck_color_map(self):

    cmap = ListedColormap(np.loadtxt("Planck_FreqMap_RGB.txt")/255.)
    cmap.set_bad("darkgray")
    cmap.set_under("white")

    return cmap

def restore_save(savfile):

    n = readsav(savfile)
    key = n.keys()
    if (len(key) != 1):
        exit(".sav file is not a combined end2end file")

    num_bands = int(n[key[0]]['num_bands'][0])
    bands = n[key[0]]['bands'][0]
    dbs_data = n[key[0]]['dbs_data_combined'][0]
    dbs_sims = n[key[0]]['dbs_sims_combined'][0] # (nsims, nspecs, nbands)

    winminell = int(n[key[0]]['winminell'][0])
    winmaxell = int(n[key[0]]['winmaxell'][0])

    winfunc_data = n[key[0]]['winfunc_data_combined'][0]
    winfunc_sims = n[key[0]]['winfunc_sims_combined'][0]

    cov_sv    = n[key[0]]['cov_sv_combined'][0]
    cov_noise = n[key[0]]['cov_noise_combined'][0]

    d = {'num_bands':num_bands, 'bands':bands, 
         'dbs_data':dbs_data, 'dbs_sims':dbs_sims,
         'winminell':winminell, 'winmaxell':winmaxell,
         'winfunc_data':winfunc_data, 'winfunc_sims':winfunc_sims,
         'cov_sv':cov_sv, 'cov_noise':cov_noise}

    return d


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

        lines1 = hp.projplot(theta, phi, 'white')
        
        return lines1


class create_sptxhfi_bandpower(object):

    def __init__(self, sptxspt_endfile=None, sptxhfi_endfile=None, hfixhfi_endfile=None, pdf_file=None, wfunc_corr=True):

        if (sptxspt_endfile is None):
            self.sptxspt_endfile = 'data/end_combined_spt150sn_spt150sn.sav'
        else:
            self.sptxspt_endfile = sptxspt_endfile

        if (sptxhfi_endfile is None):
            self.sptxhfi_endfile = 'data/end_combined_spt150sn_hfi143sn.sav'
        else:
            self.sptxhfi_endfile = sptxhfi_endfile

        if (hfixhfi_endfile is None):
            self.hfixhfi_endfile = 'data/end_combined_hfi143sn_hfi143sn.sav'
        else
            self.hfixhfi_endfile = hfixhfi_endfile

        self.wfunc_corr = wfunc_corr
        self.pdf_file = pdf_file

    def read_endfile(self):
        self.sptxspt = restore_save(self.sptxspt_endfile)
        self.sptxhfi = restore_save(self.sptxhfi_endfile)
        self.hfixhfi = restore_save(self.hfixhfi_endfile)

    def process_bandpower(self):
        ellmin = 650
        ellmax = 2500

        ip_sptxspt = np.where( (self.sptxspt['bands'] > ellmin) & (self.sptxspt['bands'] < ellmax) )[0]
        ip_sptxhfi = np.where( (self.sptxhfi['bands'] > ellmin) & (self.sptxhfi['bands'] < ellmax) )[0]
        ip_hfixhfi = np.where( (self.hfixhfi['bands'] > ellmin) & (self.hfixhfi['bands'] < ellmax) )[0]

        dbs_ave_sptxspt     = np.mean(self.sptxspt['dbs_sims'][:,1,:], axis=0)
        dbs_ave_hfixhfi     = np.mean(self.hfixhfi['dbs_sims'][:,1,:], axis=0)
        dbs_ave_sptxhfi     = np.mean(self.sptxhfi['dbs_sims'][:,1,:], axis=0)

        self.dbs_err_sptxspt     = np.sqrt(np.diag(self.sptxspt['cov_sv'][1,:,1,:]))[ip_sptxspt]
        self.dbs_err_hfixhfi     = np.sqrt(np.diag(self.hfixhfi['cov_sv'][1,:,1,:]))[ip_hfixhfi]
        self.dbs_err_sptxhfi     = np.sqrt(np.diag(self.sptxhfi['cov_sv'][1,:,1,:]))[ip_sptxhfi]

        self.dbs_data_sptxspt    = self.sptxspt['dbs_data'][1,ip_sptxspt]
        self.dbs_data_sptxhfi    = self.sptxhfi['dbs_data'][1,ip_sptxhfi]
        self.dbs_data_hfixhfi    = self.hfixhfi['dbs_data'][1,ip_hfixhfi]

        if (self.wfunc_corr):
            dbs_data_sptxhfi -= (dbs_ave_sptxhfi[ip_sptxhfi] - dbs_ave_sptxspt[ip_sptxspt])
            dbs_data_hfixhfi -= (dbs_ave_hfixhfi[ip_hfixhfi] - dbs_ave_sptxspt[ip_sptxspt])

        self.bands = self.sptxspt['bands'][ip_sptxspt]

    def plot_bandpower(self):

        fig, ax = plt.subplots()
        ax.set_position([0.1,0.1,0.85,0.75])

        ax.errorbar(self.bands, self.dbs_data_sptxspt, yerr=self.dbs_err_sptxspt, fmt='o', markersize='0', elinewidth=1.5, capsize=1.5, capthick=1.5, label=r'$\mathrm{SPT^{150}_{half1}\times\;SPT^{150}_{half2}}$')

        ax.errorbar(self.bands-12, self.dbs_data_sptxhfi, yerr=self.dbs_err_sptxhfi, fmt='o', markersize='0', elinewidth=1.5, capsize=1.5, capthick=1.5, label=r'$\mathrm{SPT^{150}_{full}\times\;HFI^{143}_{full}}$')

        ax.errorbar(self.bands+12, self.dbs_data_hfixhfi, yerr=self.dbs_err_hfixhfi, fmt='o', markersize='0', elinewidth=1.5, capsize=1.5, capthick=1.5, label=r'$\mathrm{HFI^{143}_{half1}\times\;HFI^{143}_{half2}}$')
    
        ax.legend(fontsize=12)

        if self.pdf_file is None:
            pdf_file = 'bandpower.pdf'

        plt.xlim([625,2500])
        plt.ylim([80,3000])
        plt.yscale('log')
        plt.xlabel(r'$\ell$', fontsize=16)
        plt.ylabel(r'$\mathcal{D}_{\ell}\ [\mathrm{\mu K^2}]$', fontsize=16)
        plt.savefig(pdf_file, format='pdf')
