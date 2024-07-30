import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import time
    import snpy
    import glob
    import time as systime
    import scripts.general as gen
    import numpy as np
    import matplotlib.pyplot as plt
CONSTANTS = gen.get_constants()

def snpy_csp_process():
    # paths = []
    # for path in glob.glob(CONSTANTS['csp_data_loc']+'*.txt'):
    #     if path.split('/')[-1][2:-9] in np.genfromtxt(CONSTANTS['csp_91bg_txt'], dtype=str, delimiter=', '):
    #         paths.append(path)
    # fit_args = {'save_loc': CONSTANTS['csp_saved_loc'], 'use_saved': False, 'snpy_plots': True, 'save_plots': True, 'quiet': False}
    # cspParams = gen.snpy_fit(paths=paths, **fit_args)
    cspObjs = gen.alt_data_proccesser(data_set='CSP', quiet=True)

    filter_wheel = {'B': 'B', 'H': 'H', 'J': 'J', 'Jrc2': 'Jrc2', 'V': 'V', 'Y': 'Y', 'Ydw': 'Ydw', 'g': 'g',
                    'i': 'i', 'r': 'r', 'u': 'u'}
    gen.write_ASCII(objs=cspObjs, filter_set=filter_wheel, save_loc=CONSTANTS['csp_saved_loc']+'ascii/')

    filter_wheel_colors = ['orange', 'cyan', 'violet']
    plt_args = {'color_wheel': filter_wheel_colors, 'save_loc': CONSTANTS['csp_saved_loc']+'plots/', 'y_type': 'mag', 'pause_time': 2, 'quiet': False, 'save_plots': True}
    gen.lc_plot(objs=cspObjs, **plt_args)

    # fit_args = {'save_loc': CONSTANTS['atlas_saved_loc'], 'use_saved': False, 'snpy_plots': True, 'save_plots': True, 'quiet': False}
    # atlasParams = gen.snpy_fit(paths=glob.glob(CONSTANTS['atlas_saved_loc']+'ascii/*.txt'), **fit_args)
    #
    # gen.dict_handler(choice='pack', data_dict=cspParams, path=CONSTANTS['csp_saved_loc']+'csp_saved.txt')
    #
    # gen.host_mass(CONSTANTS['csp_saved_loc']+'csp_saved.txt', save_loc='../default/', keep_data=False, update_saved=True)

    return
def snpy_atlas_process():
    atlasObjs = gen.data_proccesser(data_set='ATLAS', quiet=True)
    gen.write_ASCII(objs=atlasObjs, filter_set={'c': 'ATgr', 'o': 'ATri'}, save_loc=CONSTANTS['atlas_saved_loc']+'ascii/')

    plt_args = {'color_wheel': ['orange', 'cyan', 'violet'], 'save_loc': CONSTANTS['atlas_saved_loc']+'plots/', 'y_type': 'mag', 'pause_time': 2, 'quiet': False, 'save_plots': True}
    gen.lc_plot(objs=atlasObjs, **plt_args)

    fit_args = {'save_loc': CONSTANTS['atlas_saved_loc'], 'use_saved': False, 'snpy_plots': True, 'save_plots': True, 'quiet': False}
    # Group fit
    atlasParams = gen.snpy_fit(paths=glob.glob(CONSTANTS['atlas_saved_loc']+'ascii/*.txt'), **fit_args)
    # # Indivisual fit
    # SN = '2018ame'
    # atlasParams = gen.snpy_fit(paths=['../saved/atlas/ascii/'+SN+'_snpy.txt'], **fit_args)

    gen.dict_handler(choice='pack', data_dict=atlasParams, path=CONSTANTS['atlas_saved_loc']+'atlas_saved.txt')

    gen.host_mass(CONSTANTS['atlas_saved_loc']+'atlas_saved.txt', save_loc='../default/', keep_data=False, update_saved=True)

    return
def snpy_ztf_process():
    ztfObjs = gen.data_proccesser(data_set='ZTF', mag_unc_max=0.5, quiet=True)
    gen.write_ASCII(objs=ztfObjs, filter_set={'ZTF_g': 'g', 'ZTF_r': 'r', 'ZTF_i': 'i'}, save_loc=CONSTANTS['ztf_saved_loc']+'ascii/')

    plt_args = {'color_wheel': ['green', 'red', 'violet', 'blue'], 'save_loc': CONSTANTS['ztf_saved_loc']+'plots/', 'y_type': 'mag', 'pause_time': 2, 'quiet': False, 'save_plots': True}
    gen.lc_plot(objs=ztfObjs, **plt_args)

    fit_args = {'save_loc': CONSTANTS['ztf_saved_loc'], 'use_saved': False, 'snpy_plots': True, 'save_plots': True, 'quiet': False}
    # Group fit
    ztfParams = gen.snpy_fit(paths=glob.glob(CONSTANTS['ztf_saved_loc']+'ascii/*.txt'), **fit_args)
    # Indivisual fit
    # SN = '2018baz'
    # ztfParams = gen.snpy_fit(paths=['../saved/ztf/ascii/'+SN+'_snpy.txt'], **fit_args)

    gen.dict_handler(choice='pack', data_dict=ztfParams, path=CONSTANTS['ztf_saved_loc']+'ztf_saved.txt')

    gen.host_mass(CONSTANTS['ztf_saved_loc']+'ztf_saved.txt', save_loc='../default/', keep_data=False, update_saved=True)

    return
def snpy_indv_process(target_SN, data_set):
    if data_set == 'ATLAS':
        ascii_args = {'filter_set': {'c': 'ATgr', 'o': 'ATri'},
                      'save_loc': CONSTANTS['atlas_saved_loc']+'ascii/'}
        plt_args = {'color_wheel': ['orange', 'cyan', 'violet'], 'save_loc': CONSTANTS['atlas_saved_loc'] + 'plots/',
                    'y_type': 'mag', 'pause_time': 2, 'quiet': False, 'save_plots': True}
        fit_args = {'paths': ['../saved/ztf/ascii/'+target_SN+'_snpy.txt'], 'save_loc': CONSTANTS['ztf_saved_loc'],
                    'use_saved': False, 'snpy_plots': True, 'save_plots': True, 'quiet': False}
        save_args = {'path': CONSTANTS['ztf_saved_loc']+'ztf_saved_indv.txt'}
        replot_args = {'save_plot': False, 'stacked': False, 'save_loc': '../default/',
                       'colors': {'g': 'green', 'r': 'red', 'i': 'violet'}, 'spread': [30, 30]}
    elif data_set == 'ZTF':
        ascii_args = {'filter_set': {'ZTF_g': 'g', 'ZTF_r': 'r', 'ZTF_i': 'i'},
                      'save_loc': CONSTANTS['ztf_saved_loc']+'ascii/'}
        plt_args = {'color_wheel': ['green', 'red', 'violet', 'blue'],
                    'save_loc': CONSTANTS['ztf_saved_loc'] + 'plots/', 'y_type': 'mag', 'pause_time': 1, 'quiet': False,
                    'save_plots': True}
        fit_args = {'paths': ['../saved/ztf/ascii/'+target_SN+'_snpy.txt'], 'save_loc': CONSTANTS['atlas_saved_loc'],
                    'use_saved': True, 'snpy_plots': True, 'save_plots': True, 'quiet': False}
        save_args = {'path': CONSTANTS['ztf_saved_loc']+'ztf_saved_indv.txt'}
        replot_args = {'save_plot': False, 'stacked': False, 'save_loc': '../default/',
                       'colors': {'g': 'green', 'r': 'red', 'i': 'violet'}, 'spread': [30, 30]}
    else:
        raise ValueError("Data set must be either 'ATLAS' or 'ZTF'")
        return
    objs = gen.data_proccesser(data_set=data_set, mag_unc_max=0.5, quiet=True)
    gen.write_ASCII(objs={target_SN: objs[target_SN]}, **ascii_args)
    gen.lc_plot(objs={target_SN: objs[target_SN]} , **plt_args)
    params = gen.snpy_fit(**fit_args)
    gen.dict_handler(choice='pack', data_dict=params, **save_args)
    gen.lc_replot(CONSTANTS['ztf_saved_loc']+'models/'+target_SN+'_EBV_model2.snpy', **replot_args)
    return
def snpy_combined_process():
    # # CSP Params
    # paths = []
    # for path in glob.glob(CONSTANTS['csp_data_loc']+'*.txt'):
    #     if path.split('/')[-1][2:-9] in np.genfromtxt(CONSTANTS['csp_91bg_txt'], dtype=str, delimiter=', '):
    #         paths.append(path)
    # fit_args = {'save_loc': CONSTANTS['csp_saved_loc'], 'use_saved': True, 'snpy_plots': False, 'save_plots': False, 'quiet': False}
    # cspParams = gen.snpy_fit(paths=paths, **fit_args)
    #
    # # ATLAS-ZTF Params
    # atlasObjs = gen.data_proccesser(data_set='ATLAS', mag_unc_max=0.5, quiet=False)
    # ztfObjs = gen.data_proccesser(data_set='ZTF', mag_unc_max=0.5, quiet=False)
    # combinedObjs = atlasObjs.copy()
    # for obj in ztfObjs:
    #     if obj in combinedObjs:
    #         for cat in ['time', 'filters', 'flux', 'dflux', 'mag', 'dmag']:
    #             if cat == 'time':
    #                 combinedObjs[obj][cat] = combinedObjs[obj][cat] + 2400000 # +2400000 for MJD to JD
    #             combinedObjs[obj][cat] = np.hstack((combinedObjs[obj][cat], ztfObjs[obj][cat]))
    #     else:
    #         combinedObjs.update({obj: ztfObjs[obj]})
    # gen.write_ASCII(objs=combinedObjs, save_loc=CONSTANTS['combined_saved_loc']+'ascii/',
    #                 filter_set={'c': 'ATgr', 'o': 'ATri', 'ZTF_g': 'g', 'ZTF_r': 'r', 'ZTF_i': 'i'})
    # fit_args = {'save_loc': CONSTANTS['combined_saved_loc'], 'use_saved': True, 'snpy_plots': True, 'save_plots': False, 'quiet': False}
    # combinedParams = gen.snpy_fit(paths=glob.glob(CONSTANTS['combined_saved_loc']+'ascii/*.txt'), **fit_args)
    #
    # # Combine CSP & ATLAS-ZTF Params
    # for obj in cspParams: # Can do like this since there isn't ATLAS-ZTF & CSP overlap
    #     combinedParams.update({obj: cspParams[obj]})
    #
    # # Adding origin identifier to has determine later
    # for obj in combinedParams: # Can do like this since there isn't ATLAS-ZTF & CSP overlap
    #     if obj in cspParams:
    #         combinedParams[obj].update({'origin': 'CSP'})
    #     elif obj in atlasObjs:
    #         combinedParams[obj].update({'origin': 'ATLAS'})
    #     elif obj in ztfObjs:
    #         combinedParams[obj].update({'origin': 'ZTF'})
    #     else:
    #         print('[!!! !!!] No matches')
    #
    # gen.dict_handler(choice='pack', data_dict=combinedParams, path=CONSTANTS['combined_saved_loc']+'combined_saved.txt')

    hm_args = {'save_loc': '../default/', 'keep_data': False, 'update_saved': True, 'use_mass_key': True}
    gen.host_mass(CONSTANTS['combined_saved_loc']+'combined_saved.txt', **hm_args) # Takes a looooooooong time

    return
# ==================================================================================================================== #
def salt3_atlas_process():
    # fit_args = {'plot_data': False, 'save_plot': False}
    # atlasParams = gen.salt3_fit(data_set='ATLAS', **fit_args)
    # gen.dict_handler(data_dict=atlasParams, choice='pack', path=CONSTANTS['salt_atlas_loc']+'salt_atlas_saved.txt')
    #
    # gen.host_mass(CONSTANTS['salt_atlas_loc']+'salt_atlas_saved.txt', save_loc='../default/', keep_data=False, update_saved=True)

    # res_plot_args = {'save_loc': CONSTANTS['salt_atlas_loc'], 'sigma': [1, 1],
    #                  'labels': False, 'raw': False, 'extra_info': True, 'ignore_type': [], 'save_plot': True}
    # gen.residual_plotter(path=CONSTANTS['salt_atlas_loc']+'salt_atlas_saved.txt', x_params=['z', 'Redshift'], **res_plot_args)
    # gen.residual_plotter(path=CONSTANTS['salt_atlas_loc']+'salt_atlas_saved.txt', x_params=['host_mass', 'Host Mass'], **res_plot_args)

    return
def salt3_ztf_process():
    fit_args = {'plot_data': False, 'save_plot': False}
    atlasParams = gen.salt3_fit(data_set='ZTF', **fit_args)
    gen.dict_handler(data_dict=atlasParams, choice='pack', path=CONSTANTS['salt_ztf_loc']+'salt_ztf_saved.txt')

    gen.host_mass(CONSTANTS['salt_ztf_loc']+'salt_ztf_saved.txt', save_loc='../default/', keep_data=False, update_saved=True)

    res_plot_args = {'save_loc': CONSTANTS['salt_ztf_loc'], 'sigma': [1, 1],
                     'labels': False, 'raw': False, 'extra_info': True, 'ignore_type': [], 'save_plot': True}
    gen.residual_plotter(path=CONSTANTS['salt_ztf_loc']+'salt_ztf_saved.txt', x_params=['z', 'Redshift'], **res_plot_args)
    gen.residual_plotter(path=CONSTANTS['salt_ztf_loc']+'salt_ztf_saved.txt', x_params=['host_mass', 'Host Mass'], **res_plot_args)

    return
# ==================================================================================================================== #
def fit_reviewer():
    saved_models_loc = [glob.glob(CONSTANTS['csp_saved_loc']+'models/*.snpy'),
                        glob.glob(CONSTANTS['combined_saved_loc']+'models/*.snpy')]
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
            print('[', models.index(model)+1, '/', len(models), '] - ', objname)

            s = snpy.get_sn(model)
            gen.lc_replot(model, spread=[30, 30], stacked=False)

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
    with open(CONSTANTS['reviewed_fits_txt'], 'w') as f:
        for fits in all_fits:
            f.write(fits)
            for obj in all_fits[fits]:
                f.write(', '+"'"+obj+"'")
            f.write('\n')
    return
def plotter(choices):
    for data_set in choices:
        print("Plotting Hubble Residuals for '"+data_set+"'...")
        res_plot_args = {'save_loc': CONSTANTS[data_set+'_saved_loc'], 'sigma': [1, 1],
                         'labels': True, 'raw': False, 'extra_info': True, 'ignore_type': ['bad', 'okay'], 'save_plot': True}
        # gen.residual_plotter(path=CONSTANTS[data_set+'_saved_loc']+data_set+'_saved.txt', x_params=['z', 'Redshift'], **res_plot_args)
        # gen.residual_plotter(path=CONSTANTS[data_set+'_saved_loc']+data_set+'_saved.txt', x_params=['host_mass', 'Host Mass'], **res_plot_args)

        print("Plotting parameter histograms for '"+data_set+"'...")
        hist_args = {'save_plot': True, 'hist_bins': [100, 100, 40, 200, 80]} # ['mu', 'st', 'Tmax', 'EBVhost', 'Host Mass']
        gen.histogram_plotter(data_set, save_loc=CONSTANTS[data_set+'_saved_loc'], **hist_args)
    return

if __name__ == '__main__':
    start = systime.time()  # Runtime tracker

    # snpy_csp_process()
    # snpy_atlas_process()
    # snpy_ztf_process()
    # snpy_indv_process('2018baz', 'ZTF')
    # snpy_combined_process()

# -------------------------------------------------------------------------------------------------------------------- #

    # salt3_atlas_process()
    # salt3_ztf_process()

# -------------------------------------------------------------------------------------------------------------------- #

    # fit_reviewer()
    plotter(['combined'])

# -------------------------------------------------------------------------------------------------------------------- #

    # gen.salt3_fit(data_set='ATLAS', plot_data=False, save_plot=False)
    # gen.host_mass('../default/combined_saved.txt', save_loc='../default/', keep_data=False, update_saved=False, use_mass_key=True)

    # csp_objs = gen.alt_data_proccesser(data_set='CSP', quiet=False)
    # atlas_objs = gen.alt_data_proccesser(data_set='ATLAS', quiet=False)
    # ztf_objs = gen.alt_data_proccesser(data_set='ZTF', quiet=False)


    # # Review of good fits
    # with open(CONSTANTS['reviewed_fits_txt'], 'r') as f:
    #     paths = f.readline().split(', ')[1:]
    #     new_paths = []
    #     for p in paths:
    #         new_paths.append(p.replace("'", ''))
    # for p in new_paths:
    #     gen.lc_replot('../saved/snpy/combined/models/'+p+'_EBV_model2.snpy', spread=[30, 30], stacked=False)

    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')








