import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import snpy
    import time as systime

    import scripts.general as gen
    import scripts.csp_process as csp
    import scripts.atlas_process as atlas
    import scripts.ztf_process as ztf
    import scripts.ghost_process as ghost

""" GLOBALS """
ROOT_PATH = '../working_data/'
PROB_CHILD_TXT = ROOT_PATH+'problem_children.txt'
TNS_KEY_TXT = ROOT_PATH+'TNS_key.txt'
# COSMO_MODEL = cosmo.FlatLambdaCDM(H0=70, Om0=0.3)

if __name__ == '__main__':
    start = systime.time()  # Runtime tracker

    # csp.burns_fitting(skip_problems=False, use_saved=False, snpy_plots=True)
    # csp.ghost_host_galaxy(BURNS_SAVE_TXT, save_loc=TEST_ROOT, keep_data=False, update_saved=True)

    # atlas.atlas_collection(quiet=False, check_data=True)
    # atlas_objs = atlas.atlas_processing(err_max=1000, n_iter=0, sleep_t=5, use_TNS=True, loc_TNS=TNS_KEY_TXT)
    # gen.write_ASCII(atlas_objs, '../snpy/atlas/ascii/', quiet=True)
    atlas.atlas_snpy_fitting(n_iter=10, skip_problems=False, use_saved=False, snpy_plots=True, save_plots=True)
    # atlas.ghost_host_galaxy(ATLAS_SAVE_TXT, save_loc=TEST_ROOT, keep_data=True, update_saved=True)

    # ztf_collection(submit=True)
    # ztf_alt_collection()
    #
    # burns_plotting([]) # Options: ['reg_hist', 'res_hist', 'zvmu']
    #
    # atlas_plotting([])
    #
    # ghost_plot_args = {'plot_size': (18, 6), 'plot_ratio': [10, 1], 'hist_bins': [40, 40, 35, 35], 'labels': False, 'raw': False}
    # ghost_plotting(['altas_residualsvz'], **ghost_plot_args)
    #
    # snpy_fit_indv('2023cvq')

    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')
