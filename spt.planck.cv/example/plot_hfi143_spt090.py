import healpy as hp
import matplotlib as mpl

if __name__ == '__main__':
    if __package__ is None:
        from os import path
        import sys
        sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) )+'/source' )

        from make_figures import *
    else:
        from source.make_figures import *

    fig, ax = plt.subplots()

    file_hfi143 = 'data/HFI_SkyMap_RING_143_2048_R2.00_full.fits'
    map_hfi143 = hp.read_map(file_hfi143)

    file_spt = 'data/healpix_spt_90_nside_2048.fits'
    map_spt = hp.read_map(file_spt)
    bad_pixels = np.where(np.abs(map_spt + 1.6375e30) < 1.0e-6)


    spt090 = Healpix_SPTField(2048, background=map_spt, offset=-0.8e-4, factor=1e6, input_coord='C', bad_pixels=bad_pixels, sub=[1,2,2], cbar=False)
    orth_spt = spt090.create_background()
    spt090.create_graticule([160,140,120,100], np.arange(0,360,45))
    spt090.create_spt_area()

    hfi143 = Healpix_SPTField(2048, background=map_hfi143, offset=-0.8e-4, factor=1e6, sub=[1,2,1], cbar=False)
    orth_hfi = hfi143.create_background()
    hfi143.create_graticule([160,140,120,100], np.arange(0,360,45))
    hfi143.create_spt_area()

    cbticks = np.array([-500.0,0.0,100.0,1000.0, 1e6])
    cbticklabels = np.array([r'$-500$',r'$0$',r'$100$',r'$1000$',r'$10^6$'])
    cbticks = np.arcsinh( (cbticks-80.0) * 0.5 ) / np.log(10.0)

    cmap = hfi143.planck_color_map()
    norm = mpl.colors.Normalize(vmin=-3.1, vmax=7)

    ax.set_position([0.21, 0.10, 0.6, 0.04])
    cbar = mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='horizontal', drawedges=False)

    cbar.set_ticks(cbticks)
    cbar.ax.axes.set_xticklabels(cbticklabels, fontsize=16)
    cbar.set_label(r'$\mathrm{\mu K}$', fontsize=16)
     

    plt.savefig('hfi143_spt090.png', format='png', dpi=600)
