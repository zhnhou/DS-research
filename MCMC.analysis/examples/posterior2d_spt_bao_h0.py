import matplotlib.pyplot as plt

if __name__ == '__main__':
    if __package__ is None:
        from os import path
        import sys
        sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) )+'/source' )
        
        from make_figures import *
    else:
        from source.make_figures import *

    s12tau_chain_path = '/Users/zhenhou/Downloads/lcdm_camb_s12tau/'

    pdf_file = 'figures/posterior_spt_bao_h0.pdf'
    
    fig = Figure_CMBBAOH0('SPT', s12tau_chain_path, pdf_file=pdf_file, NoColorBar=False)
    fig.process_chain()
    fig.create_figure()
