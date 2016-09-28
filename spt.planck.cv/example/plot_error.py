if __name__ == '__main__':
    if __package__ is None:
        from os import path
        import sys
        sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) )+'/source' )

        from make_figures import *
    else:
        from source.make_figures import *

    bandpower = create_sptxhfi_bandpower(wfunc_corr=True)
    bandpower.read_endfile()
    bandpower.read_beam_cov()
    bandpower.process_bandpower()

    end_base_file    = 'data/end_combined_spt150s_hfi143s.sav'
    end_143x143_file = 'data/end_combined_hfi143sn_hfi143sn.sav'
    end_150x143_file = 'data/end_combined_spt150sn_hfi143sn.sav'

    end_150x150_file = 'data/end_combined_spt150sn_spt150sn.sav'
    beam_cov_residual_file = 'data/beam_cov_residual.sav'

    end_base = restore_end_save(end_base_file, ellmin=650, ellmax=2500)
    end_143x143 = restore_end_save(end_143x143_file, ellmin=650, ellmax=2500)
    end_150x143 = restore_end_save(end_150x143_file, ellmin=650, ellmax=2500)
    end_150x150 = restore_end_save(end_150x150_file, ellmin=650, ellmax=2500)

    end_150x143['dbs_data'] = end_base['dbs_data']

    b0 = readsav(beam_cov_residual_file)
    b1 = {}

    b1['d21d21'] = b0['d21d21'][13:50,13:50]
    b1['d31d31'] = b0['d31d31'][13:50,13:50]
    b1['d21d31'] = b0['d21d31'][13:50,13:50]
    b1['d31d21'] = b0['d31d21'][13:50,13:50]
    
    residual = create_residual_figure(end_143x143, end_150x143, end_150x150, res_beam_cov=b1)
    residual.process_end()

    
    fig, ax = plt.subplots()
    ax.set_position([0.15,0.15,0.8,0.7])

    plt.plot(bandpower.bands, bandpower.dbs_err_sptxspt, color='red', linewidth=4.5, label=r'$\Delta \mathcal{D}_b^{150\times 150}$')
    plt.plot(bandpower.bands, residual.res_info['error_150x143_150x150'], color='green', linewidth=4.5, label=r'$\Delta (\mathcal{D}_b^{150\times 143} - \mathcal{D}_b^{150\times 150})$')
    plt.plot(bandpower.bands, residual.res_info['error_143x143_150x150'], color='blue', linewidth=4.5, label=r'$\Delta (\mathcal{D}_b^{143\times 143} - \mathcal{D}_b^{150\times 150})$')
    

    ax.set_yscale('log')

    ax.set_xlim([650,2500])
    ax.set_ylim([5,400])

    ax.set_xticks([1000,1500,2000,2500])
    ax.set_xticklabels([r'$1000$',r'$1500$',r'$2000$',r'$2500$'], fontsize=22)

    ax.set_yticks([10,100])
    ax.set_yticklabels([r'$10$',r'$100$'], fontsize=22)

    plt.legend(frameon=False, fontsize=22, loc=2, labelspacing=0.18)

    ax.set_xlabel(r'$\ell$', fontsize=20)
    ax.set_ylabel(r'$\mathrm{\mu K^2}$', fontsize=20)

    plt.savefig('band_error.pdf', format='pdf', transparent=True)
