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

def restore_end_save(savfile, ellmin=None, ellmax=None):

    n = readsav(savfile)
    key = n.keys()
    if (len(key) != 1):
        exit(".sav file is not a combined end2end file")

    bands = n[key[0]]['bands'][0]

    if (not ellmin):
        ellmin = bands[0]

    if (not ellmax):
        ellmax = bands[-1]

    index1 = bands >= ellmin
    index2 = bands <= ellmax

    index = index1 * index2

    bands = bands[index]
    num_bands = len(bands)

    dbs_data = n[key[0]]['dbs_data_combined'][0][:,index]
    dbs_sims = n[key[0]]['dbs_sims_combined'][0][:,:,index] # (nsims, nspecs, nbands)

    winminell = int(n[key[0]]['winminell'][0])
    winmaxell = int(n[key[0]]['winmaxell'][0])

    winfunc_data = n[key[0]]['winfunc_data_combined'][0][:,index,:] # (nspecs, nbands, nell_wf)
    winfunc_sims = n[key[0]]['winfunc_sims_combined'][0][:,index,:]

    cov_sv    = n[key[0]]['cov_sv_combined'][0][:,index,:,index]
    cov_noise = n[key[0]]['cov_noise_combined'][0][:,index,:,index]

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

    def __init__(self, sptxspt_endfile=None, sptxhfi_endfile=None, hfixhfi_endfile=None, pdf_file=None, wfunc_corr=True, recalib=1.00):

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
        else:
            self.hfixhfi_endfile = hfixhfi_endfile

        self.wfunc_corr = wfunc_corr
        self.recalib = recalib

        if pdf_file is None:
            self.pdf_file = 'bandpower.pdf'
        else:
            self.pdf_file = pdf_file

    def read_endfile(self):
        self.sptxspt = restore_save(self.sptxspt_endfile)
        self.sptxhfi = restore_save(self.sptxhfi_endfile)
        self.hfixhfi = restore_save(self.hfixhfi_endfile)

    def read_beam_cov(self):
        tmp = readsav('data/beam_cov_150x150_150x143_143x143.sav')
        self.beam_cov_150x150 = tmp['c11'][13:50,13:50]
        self.beam_cov_150x143 = tmp['c22'][13:50,13:50]
        self.beam_cov_143x143 = tmp['c33'][13:50,13:50]

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

        self.dbs_err_sptxspt += np.sqrt(np.diag(self.beam_cov_150x150))
        self.dbs_err_sptxhfi += np.sqrt(np.diag(self.beam_cov_150x143))
        self.dbs_err_hfixhfi += np.sqrt(np.diag(self.beam_cov_143x143))

        self.dbs_data_sptxspt    = self.sptxspt['dbs_data'][1,ip_sptxspt]
        self.dbs_data_sptxhfi    = self.sptxhfi['dbs_data'][1,ip_sptxhfi]
        self.dbs_data_hfixhfi    = self.hfixhfi['dbs_data'][1,ip_hfixhfi]

        if (self.wfunc_corr):
            self.dbs_data_sptxhfi -= (dbs_ave_sptxhfi[ip_sptxhfi] - dbs_ave_sptxspt[ip_sptxspt]) * self.recalib**2
            self.dbs_data_hfixhfi -= (dbs_ave_hfixhfi[ip_hfixhfi] - dbs_ave_sptxspt[ip_sptxspt]) * self.recalib

        self.bands = self.sptxspt['bands'][ip_sptxspt]

    def plot_bandpower(self, set_yticklabels=True, set_legend=True):

        fig, ax = plt.subplots()
        ax.set_position([0.15,0.15,0.8,0.7])

        ax.errorbar(self.bands, self.dbs_data_sptxspt, yerr=self.dbs_err_sptxspt, fmt='o', markersize='0', elinewidth=1.5, capsize=1.5, capthick=1.5, label=r'$\mathrm{SPT^{150}_{half1}\times\;SPT^{150}_{half2}}$')

        ax.errorbar(self.bands-12, self.dbs_data_sptxhfi, yerr=self.dbs_err_sptxhfi, fmt='o', markersize='0', elinewidth=1.5, capsize=1.5, capthick=1.5, label=r'$\mathrm{SPT^{150}_{full}\times\;HFI^{143}_{full}}$')

        ax.errorbar(self.bands+12, self.dbs_data_hfixhfi, yerr=self.dbs_err_hfixhfi, fmt='o', markersize='0', elinewidth=1.5, capsize=1.5, capthick=1.5, label=r'$\mathrm{HFI^{143}_{half1}\times\;HFI^{143}_{half2}}$')
   
        if set_legend: 
            ax.legend(fontsize=16, frameon=False)

        plt.xlim([625,2500])
        plt.ylim([80,3000])
        plt.yscale('log')

        ax.set_xticks([1000,1500,2000,2500])
        ax.set_yticks([100,1000])

        ax.set_xticklabels([r'$1000$',r'$1500$',r'$2000$',r'$2500$'], fontsize=22)
        if set_yticklabels:
            ax.set_yticklabels([r'$10^2$', r'$10^3$'], fontsize=22)
            plt.ylabel(r'$\mathcal{D}_{\ell}\ [\mathrm{\mu K^2}]$', fontsize=22)
        else:
            ax.set_yticklabels([' ', ' '], fontsize=0)
            plt.ylabel(' ', fontsize=0)

        plt.xlabel(r'$\ell$', fontsize=22)
        plt.savefig(self.pdf_file, format='pdf', transparent=True)


class create_residual_figure(object):

    def __init__(self, end_143x143, end_150x143, end_150x150, res_beam_cov=None):

        self.end_143x143 = end_143x143
        self.end_150x143 = end_150x143
        self.end_150x150 = end_150x150


        if (res_beam_cov is None):
            self.rescale_143x143 = 1.00
            self.rescale_150x143 = 1.0114900
            self.rescale_150x150 = 1.0110200**2
        else:
            self.rescale_143x143 = 1.00
            self.rescale_150x143 = 1.0090700
            self.rescale_150x150 = 1.0090700**2

        self.res_beam_cov = res_beam_cov

        if (not np.array_equal(end_150x143['bands'], end_150x150['bands']) ):
            print "The band definition of two end2end files are different"
            exit()

        if (not np.array_equal(end_143x143['bands'], end_150x150['bands']) ):
            print "The band definition of two end2end files are different"
            exit()

    def process_end(self):
        
        end_143x143_sims_mean = self.end_143x143['dbs_sims'][:,1,:].mean(axis=0)
        end_150x143_sims_mean = self.end_150x143['dbs_sims'][:,1,:].mean(axis=0)
        end_150x150_sims_mean = self.end_150x150['dbs_sims'][:,1,:].mean(axis=0)

        winfunc_corr_150x143_150x150 = end_150x143_sims_mean - end_150x150_sims_mean
        winfunc_corr_143x143_150x150 = end_143x143_sims_mean - end_150x150_sims_mean

        res_150x143_150x150 = self.end_150x143['dbs_data'][1,:]*self.rescale_150x143 - self.end_150x150['dbs_data'][1,:]*self.rescale_150x150 - winfunc_corr_150x143_150x150
        res_143x143_150x150 = self.end_143x143['dbs_data'][1,:]*self.rescale_143x143 - self.end_150x150['dbs_data'][1,:]*self.rescale_150x150 - winfunc_corr_143x143_150x150


        res_sims_150x143_150x150 = self.end_150x143['dbs_sims'][:,1,:] - self.end_150x150['dbs_sims'][:,1,:]
        res_sims_143x143_150x150 = self.end_143x143['dbs_sims'][:,1,:] - self.end_150x150['dbs_sims'][:,1,:]

        cov_150x143_150x150 = np.cov(res_sims_150x143_150x150.transpose())
        cov_143x143_150x150 = np.cov(res_sims_143x143_150x150.transpose())

        error_150x143_150x150_nobeam = np.sqrt(np.diag(cov_150x143_150x150))
        error_143x143_150x150_nobeam = np.sqrt(np.diag(cov_143x143_150x150))

        self.error_150x143_150x150_nobeam = error_150x143_150x150_nobeam
        self.error_143x143_150x150_nobeam = error_143x143_150x150_nobeam

        if (not (self.res_beam_cov is None)):
            cov_150x143_150x150 += self.res_beam_cov['d21d21']
            cov_143x143_150x150 += self.res_beam_cov['d31d31']

        error_150x143_150x150 = np.sqrt(np.diag(cov_150x143_150x150))
        error_143x143_150x150 = np.sqrt(np.diag(cov_143x143_150x150))

        d = {'res_data_150x143_150x150':res_150x143_150x150, 'res_cov_150x143_150x150':cov_150x143_150x150,
             'res_data_143x143_150x150':res_143x143_150x150, 'res_cov_143x143_150x150':cov_143x143_150x150,
             'error_150x143_150x150':error_150x143_150x150,  'error_143x143_150x150':error_143x143_150x150}

        self.res_info = d

    def make_residual_figure(self, pdf_file):

        error_150x143_150x150 = np.sqrt(np.diag(self.res_info['res_cov_150x143_150x150']))
        error_143x143_150x150 = np.sqrt(np.diag(self.res_info['res_cov_143x143_150x150']))

        yticks = [-40, -20, 0, 20, 40]
        xticks = [650,1000,1500,2000,2500]

        fig = plt.figure()
        ax1 = fig.add_subplot(211)
        ax1.set_position([0.13,0.50,0.85,0.25])

        ax1.plot([0,3000],[0,0], color='black', linewidth=0.5, zorder=0)
        ax1.errorbar(self.end_150x143['bands'], self.res_info['res_data_150x143_150x150'], yerr=error_150x143_150x150, fmt='o', markersize='0', elinewidth=2, capsize=2., capthick=2., zorder=3)
#ax1.errorbar(self.end_150x143['bands'], self.res_info['res_data_150x143_150x150'], yerr=self.error_150x143_150x150_nobeam, fmt='o', markersize='0', elinewidth=2, capsize=0., capthick=2.)
        
        ax1.set_xlim([625,2525])
        ax1.set_ylim([-55,55])
        ax1.set_xticks(xticks)
        ax1.set_yticks(yticks)
        ax1.axes.set_xticklabels([" "," "," "," "," "])
        ax1.axes.set_yticklabels(["$-40$","$-20$","$0$","$20$","$40$"], fontsize=16)

        ax1.text(750, 35, r"$\mathcal{D}_b^{150 \times 143} - \mathcal{D}_b^{150 \times 150}$", fontsize=18)


        yticks = [-100, -50, 0, 50, 100]
        ax2 = fig.add_subplot(212)
        ax2.set_position([0.13,0.25,0.85,0.25])
        
        ax2.plot([0,3000],[0,0], color='black', linewidth=0.5, zorder=0)
        ax2.errorbar(self.end_143x143['bands'], self.res_info['res_data_143x143_150x150'], yerr=error_143x143_150x150, fmt='o', markersize='0', elinewidth=2., capsize=2., capthick=2., zorder=3)
        ax2.set_xlim([625,2525])
        ax2.set_ylim([-137.5,137.5])

        ax2.set_xticks(xticks)
        ax2.set_yticks(yticks)
        ax2.axes.set_xticklabels(["$650$","$1000$","$1500$","$2000$","$2500$"], fontsize=16)
        ax2.axes.set_yticklabels(["$-100$","$-50$","$0$","$50$","$100$"], fontsize=16)
        ax2.set_xlabel("$\ell$", fontsize=20)

        ax2.text(750, 87.5, r"$\mathcal{D}_b^{143 \times 143} - \mathcal{D}_b^{150 \times150}$", fontsize=18)
        ax2.text(400,137.5,"$\Delta \mathcal{D}_b\,[\mathrm{\mu K^2}]$", rotation=90, ha='center', va='center', fontsize=20)
        
        plt.savefig(pdf_file, format='pdf')
        plt.clf()
