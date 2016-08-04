import matplotlib.pyplot as plt

if __name__ == '__main__':
    if __package__ is None:
        from os import path
        import sys
        sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) ) )
        
        from mcmc_analysis import *
    else:
        from ..mcmc_analysis import *

    s12tau_chain_path = '/Users/zhenhou/Downloads/lcdm_camb_s12tau/'
    feature = ['H0', 'rs_zdrag', 'Dv_0.57']

    m = mcmc_analysis(chain_path=s12tau_chain_path, feature=feature)
    m.chain['Dv057rs'] = 100.0 * m.chain['rs_zdrag'] / m.chain['Dv_0.57']

    post2d = m.create_2d_posterior('H0', 'Dv057rs', xra=None, yra=None, nbins_raw=80, frate=2.0)
    
    fig, ax = plt.subplots()
    contour = plt.contour(post2d['grid_x'], post2d['grid_y'], post2d['posterior2d'], post2d['levels'], linewidths=2)

    ax.set_aspect(np.float(post2d['dx_raw'] / post2d['dy_raw']))

    plt.show()
