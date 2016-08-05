import matplotlib.pyplot as plt
from scipy.special import *

def calc_clevel_gmm(logprob):
    z = np.exp(-logprob)
    nm = np.sum(z)
    z /= nm

    psort = np.sort(z.flatten())[::-1]
    pcum  = np.cumsum(psort)

    clevels = []
    sigma = [2.0, 1.0]
    ptes = erf(sigma / np.sqrt(2.0))

    for p in ptes:
        ind = np.where(pcum > p)[0][0]
        clevels.append(psort[ind])

    return z, clevels


if __name__ == '__main__':
    if __package__ is None:
        from os import path
        import sys
        sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) ) )
        
        from mcmc_analysis import *
    else:
        from ..mcmc_analysis import *

    s12tau_chain_path = '/Users/zhenhou/Downloads/lcdm_camb_s12tau/'

    
    x_param = 'ns'
    x_param = 'H0'

    y_param = 'Omegabh2'
    y_param = 'Dv057rs'

    z_param = 'Omegamh2'

    feature = [x_param, 'Dv_0.57', 'rs_zdrag', z_param]

    m = mcmc_analysis(chain_path=s12tau_chain_path, feature=feature)
    m.chain['Dv057rs'] = 100.0 * m.chain['rs_zdrag'] / m.chain['Dv_0.57']
#    m.chain['ns'] *= 100.0
#    m.chain['Omegabh2'] *= 1000.0

    x_select, y_select, z_select = m.create_2d_scatter(x_param, y_param, z_param, num_select=10000)
    scatter = plt.scatter(x_select, y_select, c=z_select, s=1, vmin=0.105, vmax=0.155, edgecolors='face', cmap=plt.get_cmap('YlGnBu_r'))

    X, Y, Z = m.create_2d_gmm(x_param, y_param)

    z, clevels = calc_clevel_gmm(Z)

    plt.contour(X, Y, z, levels=clevels, colors='k')
    plt.axes().set_aspect(11)
    plt.show()
