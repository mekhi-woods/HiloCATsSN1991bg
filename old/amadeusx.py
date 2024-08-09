import time
import warnings

import numpy as np
import matplotlib.pyplot as plt

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

def depricated():
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
    # # Combined ATLAS-ZTF
    # import numpy as np
    # atlasObjs = gen.data_proccesser(data_set='ATLAS', quiet=True)
    # ztfObjs = gen.data_proccesser(data_set='ZTF', mag_unc_max=0.5, quiet=True)
    # combinedObjs = atlasObjs.copy()
    # for obj in ztfObjs:
    #     if obj in combinedObjs:
    #         for cat in ['time', 'filters', 'flux', 'dflux', 'mag', 'dmag']:
    #             if cat == 'time':
    #                 combinedObjs[obj][cat] = combinedObjs[obj][cat] + 2400000 # +2400000 for MJD to JD
    #             combinedObjs[obj][cat] = np.hstack((combinedObjs[obj][cat], ztfObjs[obj][cat]))
    #     else:
    #         combinedObjs.update({obj: ztfObjs[obj]})
    #
    # plt_args = {'color_wheel': ['green', 'red', 'violet', 'blue', 'orange', 'yellow'],
    #             'save_loc': '../snpy/combined/plots/', 'y_type': 'mag', 'pause_time': 1, 'quiet': False, 'save_plots': True}
    # gen.lc_plot(objs=combinedObjs, **plt_args)
    #
    # gen.write_ASCII(objs=combinedObjs, save_loc='../snpy/combined/ascii/',
    #                 filter_set={'c': 'ATgr', 'o': 'ATri', 'ZTF_g': 'g', 'ZTF_r': 'r', 'ZTF_i': 'i'})
    #
    # fit_args = {'save_loc': '../snpy/combined/', 'use_saved': False, 'snpy_plots': True, 'save_plots': False, 'quiet': False}
    # combinedParams = gen.snpy_fit(paths=glob.glob('../snpy/combined/ascii/*.txt'), **fit_args)
    # gen.dict_packer(combinedParams, '../snpy/combined/combined_saved.txt')
    #
    # ghost.ghost_host_galaxy('../snpy/combined/combined_saved.txt', save_loc='../tests/', keep_data=False, update_saved=True)
    #
    # bad_fits = ['']
    # res_plot_args = {'labels': True, 'raw': False, 'extra_info': True, 'bad_fits': bad_fits}
    # gen.residual_ploter(path='../snpy/combined/combined_saved.txt', x_params=['z', 'Redshift'], **res_plot_args)
    # gen.residual_ploter(path='../snpy/combined/combined_saved.txt', x_params=['logstellarmass', 'Host Mass'], labels=False)

    return
def csp_process():
    paths = []
    for path in glob.glob('../data/CSPdata/*.txt'):
        if path.split('/')[-1][2:-9] in np.genfromtxt('../txts/91bglike_justnames.txt', dtype=str, delimiter=', '):
            paths.append(path)
    fit_args = {'save_loc': '../snpy/csp/', 'use_saved': False, 'snpy_plots': True, 'save_plots': True, 'quiet': False}
    cspParams = gen.snpy_fit(paths=paths, **fit_args)
    return
def atlas_process():
    atlasObjs = gen.data_proccesser(data_set='ATLAS', quiet=False)
    gen.write_ASCII(objs=atlasObjs, filter_set={'c': 'ATgr', 'o': 'ATri'}, save_loc='../snpy/atlas/ascii/')

    # plt_args = {'color_wheel': ['orange', 'cyan', 'violet'], 'save_loc': '../snpy/atlas/plots/', 'y_type': 'mag', 'pause_time': 2, 'quiet': False, 'save_plots': True}
    # gen.lc_plot(objs=atlasObjs, **plt_args)

    fit_args = {'save_loc': '../snpy/atlas/', 'use_saved': False, 'snpy_plots': True, 'save_plots': False, 'quiet': False}
    atlasParams = gen.snpy_fit(paths=glob.glob('../snpy/atlas/ascii/*.txt'), **fit_args)
    gen.dict_handler(choice='pack', data_dict=atlasParams, '../snpy/atlas/atlas_saved.txt')
    return
def ztf_process():
    ztfObjs = gen.data_proccesser(data_set='ZTF', mag_unc_max=0.5, quiet=True)
    gen.write_ASCII(objs=ztfObjs, filter_set={'ZTF_g': 'g', 'ZTF_r': 'r', 'ZTF_i': 'i'}, save_loc='../snpy/ztf/ascii/')

    # plt_args = {'color_wheel': ['green', 'red', 'violet', 'blue'], 'save_loc': '../snpy/ztf/plots/', 'y_type': 'mag', 'pause_time': 1, 'quiet': False, 'save_plots': True}
    # gen.lc_plot(objs=ztfObjs['2018lph'], **plt_args)

    fit_args = {'save_loc': '../snpy/ztf/', 'use_saved': False, 'snpy_plots': True, 'save_plots': True, 'quiet': False}
    ztfParams = gen.snpy_fit(paths=glob.glob('../snpy/ztf/ascii/*.txt'), **fit_args)
    gen.dict_handler(choice='pack', data_dict=ztfParams, path='../snpy/ztf/ztf_saved.txt')
    return
def combined_process():
    # CSP Params
    paths = []
    for path in glob.glob('../data/CSPdata/*.txt'):
        if path.split('/')[-1][2:-9] in np.genfromtxt('../txts/91bglike_justnames.txt', dtype=str, delimiter=', '):
            paths.append(path)
    fit_args = {'save_loc': '../snpy/csp/', 'use_saved': False, 'snpy_plots': True, 'save_plots': True, 'quiet': False}
    cspParams = gen.snpy_fit(paths=paths, **fit_args)

    # ATLAS-ZTF Params
    atlasObjs = gen.data_proccesser(data_set='ATLAS', mag_unc_max=0.5, quiet=False)
    ztfObjs = gen.data_proccesser(data_set='ZTF', mag_unc_max=0.5, quiet=False)
    combinedObjs = atlasObjs.copy()
    for obj in ztfObjs:
        if obj in combinedObjs:
            for cat in ['time', 'filters', 'flux', 'dflux', 'mag', 'dmag']:
                if cat == 'time':
                    combinedObjs[obj][cat] = combinedObjs[obj][cat] + 2400000 # +2400000 for MJD to JD
                combinedObjs[obj][cat] = np.hstack((combinedObjs[obj][cat], ztfObjs[obj][cat]))
        else:
            combinedObjs.update({obj: ztfObjs[obj]})
    gen.write_ASCII(objs=combinedObjs, save_loc='../snpy/combined/ascii/',
                    filter_set={'c': 'ATgr', 'o': 'ATri', 'ZTF_g': 'g', 'ZTF_r': 'r', 'ZTF_i': 'i'})
    fit_args = {'save_loc': '../snpy/combined/', 'use_saved': False, 'snpy_plots': True, 'save_plots': True, 'quiet': False}
    combinedParams = gen.snpy_fit(paths=glob.glob('../snpy/combined/ascii/*.txt'), **fit_args)

    # Combine CSP & ATLAS-ZTF Params
    for obj in cspParams: # Can do like this since there isn't ATLAS-ZTF & CSP overlap
        combinedParams.update({obj: cspParams[obj]})

    gen.dict_handler(choice='pack', data_dict=combinedParams, path='../snpy/combined/combined_saved.txt')

    ghost.ghost_host_galaxy('../snpy/combined/combined_saved.txt', save_loc='../tests/', keep_data=False, update_saved=True) # Takes a looooooooong time

    return
def fit_reviewer():
    saved_models_loc = [glob.glob('../snpy/csp/*.snpy'),
                        glob.glob('../snpy/combined/*.snpy')]
    all_bad_fit, all_okay_fit, all_good_fit = [], [], []
    for models in saved_models_loc:
        print('========================================================================')
        print('Fits for', models[0].split('/')[-2].upper(), 'data...')
        print('========================================================================')
        bad_fit, okay_fit, good_fit = [], [], []
        for model in models:
            objname = model.split('/')[-1].split('_')[0]
            if objname[:2] == 'SN':
                objname = objname[2:]
            print('------------------------------------------------------------------------')
            print(objname)

            s = snpy.get_sn(model)
            s.plot()
            plt.show()

            selecting = True
            while selecting:
                selection = input('How is this fit? [1-bad, 2-okay, 3-good]: ')
                if selection == '1':
                    bad_fit.append(objname)
                    selecting = False
                elif selection == '2':
                    okay_fit.append(objname)
                    selecting = False
                elif selection == '3':
                    good_fit.append(objname)
                    selecting = False
                else:
                    print('[!!!] Unknown selection.')
        print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print('Good fits: ', good_fit)
        print('Okay fits: ', okay_fit)
        print('Bad fits: ', bad_fit)

        all_good_fit = all_good_fit + good_fit
        all_okay_fit = all_okay_fit + okay_fit
        all_bad_fit = all_bad_fit + bad_fit


    all_fits = {'good': all_good_fit, 'okay': all_okay_fit, 'bad': all_bad_fit}
    with open('../reviewed_fits.txt', 'w') as f:
        for fits in all_fits:
            f.write(fits)
            for obj in all_fits[fits]:
                f.write(', '+"'"+obj+"'")
            f.write('\n')
    return
def plotter():
    # Residual Plotter
    res_plot_args = {'labels': True, 'raw': False, 'extra_info': True, 'ignore_type': []}
    gen.residual_ploter(path='../snpy/ztf/ztf_saved.txt', x_params=['z', 'Redshift'], **res_plot_args)
    # gen.residual_ploter(path='../snpy/ztf/ztf_saved.txt', x_params=['logstellarmass', 'Host Mass'], **res_plot_args)

    # Histograms
    # objs = np.genfromtxt('../snpy/combined/combined_saved.txt', dtype=str, skip_header=1, delimiter=', ')
    # num_objs = len(objs[:, 0])
    # plot_titles = ['mu', 'st', 'Tmax', 'EBVhost', 'Host Mass']
    # hist_indexs = [7, 8, 9, 10, 15]
    # hist_err_indexs = [11, 12, 13, 14, 16]
    # hist_bins = [40, 40, 40, 40, 80]
    # for i in range(len(hist_indexs)):
    #     plt.figure(figsize=(12, 6))
    #     plt.hist(objs[:, hist_indexs[i]].astype(float), bins=hist_bins[i])
    #     plt.title(plot_titles[i]+'\nNumber of Transients: '+str(num_objs))
    #     plt.show()

    # import numpy
    # import pyplot
    #
    # x = numpy.random.rand(1000)
    # y, bin_edges = numpy.histogram(x, bins=10)
    # bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    #
    # pyplot.errorbar(
    #     bin_centers,
    #     y,
    #     yerr=y ** 0.5,
    #     marker='.',
    #     drawstyle='steps-mid'
    # )
    # pyplot.show()
    return
def test():
    SN = '2020mfd'

    ztfObjs = gen.data_proccesser(data_set='ZTF', mag_unc_max=0.5, quiet=True)
    gen.write_ASCII(objs={SN: ztfObjs[SN]}, filter_set={'ZTF_g': 'g', 'ZTF_r': 'r', 'ZTF_i': 'i'}, save_loc='../tests/ztf_testing/')

    plt_args = {'color_wheel': ['green', 'red', 'violet', 'blue'], 'save_loc': '../tests/ztf_testing/', 'y_type': 'mag', 'pause_time': 1, 'quiet': False, 'save_plots': True}
    gen.lc_plot(objs={SN: ztfObjs[SN]}, **plt_args)

    fit_args = {'save_loc': '../tests/ztf_testing/', 'use_saved': False, 'snpy_plots': True, 'save_plots': True, 'quiet': False}
    ztfParams = gen.snpy_fit(paths=['../tests/ztf_testing/'+SN+'_snpy.txt'], **fit_args)

    hm_args = {'save_loc': '../tests/host_mass_tings/', 'keep_data': False, 'update_saved': True}
    gen.host_mass('../tests/host_mass_tings/atlas_saved.txt', **hm_args)
    return

if __name__ == '__main__':
    start = systime.time()  # Runtime tracker

    # csp_process()
    # atlas_process()
    # ztf_process()
    # combined_process()
    # fit_reviewer()
    plotter()

    # test()

    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')








