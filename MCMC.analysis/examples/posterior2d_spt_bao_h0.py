import matplotlib.pyplot as plt
from make_figures import *



if __name__ == '__main__':
    if __package__ is None:
        from os import path
        import sys
        sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) ) )
        
        from mcmc_analysis import *
    else:
        from ..mcmc_analysis import *

    s12tau_chain_path = '/Users/zhenhou/Downloads/lcdm_camb_s12tau/'

    pdf_file = '../../figures/posterior_spt_bao_h0.pdf'
    
    x_param = 'H0'
    y_param = 'Dv057rs'
    z_param = 'Omegamh2'

    feature = [x_param, 'Dv_0.57', 'rs_zdrag', z_param]

    m = mcmc_analysis(chain_path=s12tau_chain_path, feature=feature)
    m.chain['Dv057rs'] = 100.0 * m.chain['rs_zdrag'] / m.chain['Dv_0.57']

    fp = FigProp_CMBBAOH0()

    fig, ax = plt.subplots()

    plt.rc('font', family='serif')

    x_select, y_select, z_select = m.create_2d_scatter(x_param, y_param, z_param, num_select=10000)
    scatter = plt.scatter(x_select, y_select, c=z_select, s=1, vmin=0.105, vmax=0.155, edgecolors='face', cmap=plt.get_cmap('YlGnBu_r'))

    X, Y, Z = m.create_2d_gmm(x_param, y_param)

    z, clevels = calc_clevel_gmm(Z)

    plt.contour(X, Y, z, levels=clevels, colors='k')

    ax.set_xlim([fp.xmin, fp.xmax])
    ax.set_ylim([fp.ymin, fp.ymax])

    ax.set_xticks(fp.xticks)
    ax.set_yticks(fp.yticks)
    
    ax.set_xlabel(fp.xlabel, fontsize=16)
    ax.set_ylabel(fp.ylabel, fontsize=16)

    ax.set_xticklabels(fp.xticklabels, fontsize=16)
    ax.set_yticklabels(fp.yticklabels, fontsize=16)

    ax.set_aspect((fp.xmax-fp.xmin)/(fp.ymax-fp.ymin))

    plt.savefig(pdf_file, format='pdf')
