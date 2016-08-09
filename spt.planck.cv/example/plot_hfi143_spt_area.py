import healpy as hp

if __name__ == '__main__':
    if __package__ is None:
        from os import path
        import sys
        sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) )+'/source' )

        from make_figures import *
    else:
        from source.make_figures import *

    file_hfi143 = 'data/HFI_SkyMap_RING_143_2048_R2.00_full.fits'
    map_hfi143 = hp.read_map(file_hfi143)

    spt = Healpix_SPTField(2048, background=map_hfi143, offset=-0.8e-4, factor=1e6)
    spt.create_background()

