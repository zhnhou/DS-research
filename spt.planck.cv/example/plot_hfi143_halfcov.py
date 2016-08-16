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


    hm1 = Healpix_MollSPT(m1, vmin=0, vmax=4000, sub=[2,1,1], title=' ', cbar=False, margins=[0.1,0,0.1,0])
    figure1 = hm1.create_background()
    lines1 = hm1.create_spt_area()

    hm2 = Healpix_MollSPT(m2, vmin=0, vmax=4000, sub=[2,1,2], title=' ', cbar=False, margins=[0.1,0.04,0.1,-0.04])
    figure2 = hm2.create_background()
    lines2 = hm2.create_spt_area()

    plt.show()
