import healpy as hp

if __name__ == '__main__':
    if __package__ is None:
        from os import path
        import sys
        sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) )+'/source' )

        from make_figures import *
    else:
        from source.make_figures import *

    file_spt = 'data/healpix_spt_90_nside_2048.fits'
    map_spt = hp.read_map(file_spt)

    bad_pixels = np.where(np.abs(map_spt + 1.6375e30) < 1.0e-6)

    spt = Healpix_SPTField(2048, background=map_spt, offset=-0.8e-4, factor=1e6, input_coord='C', bad_pixels=bad_pixels)
    spt.create_background()
