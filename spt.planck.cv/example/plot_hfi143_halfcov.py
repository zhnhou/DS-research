import matplotlib as mpl
import heapy as hp

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

    m1 = 
