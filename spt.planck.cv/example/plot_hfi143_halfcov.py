import matplotlib as mpl
import healpy as hp

if __name__ == '__main__':
    if __package__ is None:
        from os import path
        import sys
        sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) )+'/source' )

        from make_figures import *
    else:
        from source.make_figures import *


    file_half1 = 'data/HFI_CovMap_RING_143_2048_R2.00_halfmission-1.fits'
    file_half2 = 'data/HFI_CovMap_RING_143_2048_R2.00_halfmission-2.fits'

    m1 = hp.read_map(file_half1)*1e12
    m2 = hp.read_map(file_half2)*1e12

    fig, ax = plt.subplots()

    hm1 = Healpix_MollSPT(m1, vmin=0, vmax=4000, sub=[2,1,1], title=r'$\mathrm{Halfmission 1}$', cbar=False, margins=[0.112,0.008,0.112,-0.008])
    figure1 = hm1.create_background()
    lines1 = hm1.create_spt_area()

    hm2 = Healpix_MollSPT(m2, vmin=0, vmax=4000, sub=[2,1,2], title=r'$\mathrm{Halfmission 2}$', cbar=False, margins=[0.112,0.046,0.112,-0.046])
    figure2 = hm2.create_background()
    lines2 = hm2.create_spt_area() 

    cmap = plt.get_cmap('jet')
    norm = mpl.colors.Normalize(vmin=0, vmax=4000)

    ax.set_position([0.25, 0.1, 0.5, 0.02])
    cbar = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='horizontal', drawedges=False)

    cbticks = [0, 1000, 2000, 3000, 4000]
    cbticklabels = [r'$0$', r'$1000$', r'$2000$', r'$3000$', r'$4000$']
    cbar.set_ticks(cbticks)
    cbar.ax.axes.set_xticklabels(cbticklabels, fontsize=14)
    cbar.set_label(r'$\mathrm{\mu K}^2$', fontsize=14)
    ax.set_zorder(5)

    plt.savefig('hfi143cov_half.pdf', format='pdf', dpi=600) 
