import scripts.general as gen
import scripts.csp_process as csp

""" GLOBALS """
tns_bot_id, tns_bot_name, tns_bot_api_key = '73181', 'YSE_Bot1', '0d771345fa6b876a5bb99cd5042ab8b5ae91fc67'
ROOT_PATH = 'working_data/'
PROB_CHILD_TXT = ROOT_PATH+'problem_children.txt'
TNS_KEY_TXT = ROOT_PATH+'TNS_key.txt'

COSMO_MODEL = cosmo.FlatLambdaCDM(H0=70, Om0=0.3)


if __name__ == '__main__':
    start = systime.time()  # Runtime tracker

    gen.read_DR3()



    # csp.burns_fitting(skip_problems=False, use_saved=False, snpy_plots=False)
    # csp.ghost_host_galaxy(BURNS_SAVE_TXT, save_loc=TEST_ROOT, keep_data=False, update_saved=True)

    # atlas_collection(quiet=False, check_data=True)
    # atlas_objs = atlas_processing(err_max=1000, n_iter=0, sleep_t=5, use_TNS=True)
    # write_ASCII(atlas_objs, SNPY_ATLAS_ASCII, quiet=True)
    # atlas_snpy_fitting(n_iter=0, skip_problems=False, use_saved=False, snpy_plots=True, save_plots=True)
    # ghost_host_galaxy(ATLAS_SAVE_TXT, save_loc=TEST_ROOT, keep_data=True, update_saved=True)

    # ztf_collection(submit=True)
    # ztf_alt_collection()

    # burns_plotting([]) # Options: ['reg_hist', 'res_hist', 'zvmu']

    # atlas_plotting([])

    # ghost_plot_args = {'plot_size': (18, 6), 'plot_ratio': [10, 1], 'hist_bins': [40, 40, 35, 35], 'labels': False, 'raw': False}
    # ghost_plotting(['altas_residualsvz'], **ghost_plot_args)

    # snpy_fit_indv('2023cvq')

    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')
