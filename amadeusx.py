import time
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import snpy
    import glob
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

if __name__ == '__main__':
    start = systime.time()  # Runtime tracker

    # csp.burns_fitting(skip_problems=False, use_saved=False, snpy_plots=True)
    # ghost.ghost_host_galaxy('../snpy/burns/burns_saved.txt', save_loc='../tests/', keep_data=False, update_saved=False)

    # atlas.atlas_collection(quiet=False, check_data=True)
    # atlas_objs = atlas.atlas_processing(err_max=1000, n_iter=0, sleep_t=5, use_TNS=True, loc_TNS=TNS_KEY_TXT)
    # gen.write_ASCII(atlas_objs, '../snpy/atlas/ascii/', quiet=True)
    # atlas.atlas_snpy_fitting(n_iter=0, skip_problems=False, use_saved=False, snpy_plots=True, save_plots=True)
    # ghost.ghost_host_galaxy('../snpy/atlas/atlas_saved.txt', save_loc='../tests/', keep_data=False, update_saved=True)

    # ztf.ztf_wget(submit=True)
    # ZTFobjs = ztf.ztf_alt_processing(paths = glob.glob('../data/ZTF/*.txt'), min_pts = 10, mag_err_max = 0.75, flux_err_max = 80)
    # ztf.ztf_write_ASCII(ZTFobjs, '../snpy/ztf/ascii/', quiet=False)
    # ztf.ztf_snpy_fitting(n_iter=0, skip_problems=False, use_saved=False, snpy_plots=True, save_plots=True)

    # burns_plotting([]) # Options: ['reg_hist', 'res_hist', 'zvmu']
    # atlas.atlas_plotting(['reg_hist'])
    # ztf.ztf_plotting(ZTFobjs, choice = 'flux', zoom=0, plot_filters=['g', 'r', 'i'], sigma=1, pause_time=3)
    # ztf.ztf_plotting(ZTFobjs, choice = 'mag', zoom=0, plot_filters=['g', 'r', 'i'], sigma=0, pause_time=3)

    # ghost.combined_data()

    # ghost.atlas_ztf_joint_fitting()

    # ghost_plot_args = {'plot_size': (18, 6), 'plot_ratio': [10, 1], 'hist_bins': [40, 40, 35, 35], 'labels': True, 'raw': False, 'save': True, 'ignore_type': ['bad']}
    # ghost.ghost_plotting(['combined_residualsvmass'], **ghost_plot_args)

    # snpy_fit_indv('2023cvq')

    # ztfObjs = gen.data_proccesser(data_set='ZTF', flux_err_max=100, mag_err_max=20, quiet=True)

    # ================================================================================================================ #

    # ATLAS Process
    # atlasObjs = gen.data_proccesser(data_set='ATLAS', quiet=True)
    # gen.write_ASCII(objs=atlasObjs, data_set='ATLAS', save_loc='../snpy/atlas/ascii/')
    #
    # plt_args = {'color_wheel': ['orange', 'cyan', 'violet'], 'save_loc': '../snpy/atlas/plots/', 'y_type': 'mag', 'pause_time': 2, 'quiet': False, 'save_plots': True}
    # gen.lc_plot(objs=atlasObjs, **plt_args)

    # fit_args = {'save_loc': '../snpy/atlas/', 'use_saved': True, 'snpy_plots': False, 'save_plots': False, 'quiet': False}
    # atlasParams = gen.snpy_fit(paths=glob.glob('../snpy/atlas/ascii/*.txt'), **fit_args)
    # gen.dict_packer(atlasParams, '../snpy/atlas/atlas_saved.txt')

    # ================================================================================================================ #

    # ZTF Process
    # ztfObjs = gen.data_proccesser(data_set='ZTF', mag_unc_max=0.5, quiet=True)
    # gen.write_ASCII(objs=ztfObjs, data_set='ZTF', save_loc='../snpy/ztf/ascii/')
    #
    # plt_args = {'color_wheel': ['green', 'red', 'violet', 'blue'], 'save_loc': '../snpy/ztf/plots/', 'y_type': 'mag', 'pause_time': 1, 'quiet': False, 'save_plots': True}
    # gen.lc_plot(objs=ztfObjs, **plt_args)
    #
    # fit_args = {'save_loc': '../snpy/ztf/', 'use_saved': True, 'snpy_plots': True, 'save_plots': False, 'quiet': False}
    # atlasParams = gen.snpy_fit(paths=glob.glob('../snpy/ztf/ascii/*.txt'), **fit_args)
    # gen.dict_packer(atlasParams, '../snpy/atlas/ztf_saved.txt')

    # ================================================================================================================ #

    import numpy as np
    # Combined ATLAS-ZTF
    atlasObjs = gen.data_proccesser(data_set='ATLAS', quiet=True)
    ztfObjs = gen.data_proccesser(data_set='ZTF', mag_unc_max=0.5, quiet=True)
    combinedObjs = atlasObjs.copy()
    for obj in ztfObjs:
        if obj in combinedObjs:
            for cat in ['time', 'filters', 'flux', 'dflux', 'mag', 'dmag']:
                if cat == 'time':
                    combinedObjs[obj][cat] = combinedObjs[obj][cat] + 2400000 # +2400000 for MJD to JD
                combinedObjs[obj][cat] = np.hstack((combinedObjs[obj][cat], ztfObjs[obj][cat]))
        # else:
        #     combinedObjs.update({obj: ztfObjs[obj]})

    gen.write_ASCII(objs=combinedObjs, data_set='ZTF', save_loc='../snpy/ztf/ascii/')

    # plt_args = {'color_wheel': ['green', 'red', 'violet', 'blue', 'orange', 'yellow'],
    #             'save_loc': '../snpy/combined/plots/', 'y_type': 'mag', 'pause_time': 1, 'quiet': False, 'save_plots': True}
    # gen.lc_plot(objs=combinedObjs, **plt_args)

    # fit_args = {'save_loc': '../snpy/ztf/', 'use_saved': True, 'snpy_plots': True, 'save_plots': False, 'quiet': False}
    # atlasParams = gen.snpy_fit(paths=glob.glob('../snpy/ztf/ascii/*.txt'), **fit_args)
    # gen.dict_packer(atlasParams, '../snpy/atlas/ztf_saved.txt')



    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')








