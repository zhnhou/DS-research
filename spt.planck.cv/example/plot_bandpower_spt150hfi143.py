if __name__ == '__main__':
    if __package__ is None:
        from os import path
        import sys
        sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) )+'/source' )

        from make_figures import *
    else:
        from source.make_figures import *

    pdf_file = 'bandpower_spt150hfi143.pdf'

    bp = create_sptxhfi_bandpower(pdf_file=pdf_file, wfunc_corr=False)
    bp.read_endfile()
    bp.read_beam_cov()
    bp.process_bandpower()
    bp.plot_bandpower()

    del bp

    pdf_file = 'bandpower_spt150hfi143_wfuncCorr_recalib.pdf'

    bp = create_sptxhfi_bandpower(pdf_file=pdf_file, wfunc_corr=True, recalib=1.0088)
    bp.read_endfile()
    bp.read_beam_cov()
    bp.process_bandpower()
    bp.plot_bandpower(set_yticklabels=False, set_legend=False)

    del bp

    pdf_file = 'bandpower_spt150hfi143_wfuncCorr_recalib_large.pdf'

    bp = create_sptxhfi_bandpower(pdf_file=pdf_file, wfunc_corr=True, recalib=1.0088)
    bp.read_endfile()
    bp.read_beam_cov()
    bp.process_bandpower()
    bp.plot_bandpower_large()
