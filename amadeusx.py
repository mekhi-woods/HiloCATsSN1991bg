import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import glob
    import time as systime
    import scripts.general as gen
    import scripts.snpyfit as snpy_fit
    import scripts.saltfit as salt_fit
    import scripts.plotter as pltr

CONSTANTS = gen.get_constants()

def fit_reviewer(choice='COMBINED'):
    if choice == 'COMBINED':
        saved_models_loc = [glob.glob(CONSTANTS['csp_saved_loc'] + 'models/*.snpy'),
                            glob.glob(CONSTANTS['combined_saved_loc'] + 'models/*.snpy')]
    else:
        raise ValueError("[!!!] No matches! ['COMBINED']")

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

            pltr.lc_replot(model, spread=[30, 30], stacked=False)

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
def fit_SNooPy(choice='CSP', single=''):
    snpy_fit_params = {'mag_unc_max': 0, 'flux_unc_max': 0,
                       'replot': False, 'use_saved': True, 'show_plots': True, 'save': True, 'quiet': False}
    if len(single) > 0:
        snpy_fit.indivisual_process(SN=single, data_set=choice, **snpy_fit_params)
    elif choice in ['ATLAS', 'ZTF', 'CSP']:
        snpy_fit.dataset_process(data_set=choice, **snpy_fit_params)
    elif choice in ['COMBINED']:
        snpy_fit.combined_process(**snpy_fit_params)
    return
def plot_SNooPy(choice='CSP', types=['res', 'hist']):
    path = CONSTANTS[choice.lower()+'_saved_loc'] + choice.lower()+'_saved.txt'
    if 'res' in types:
        res_plt_args = {'sigma': [1, 1], 'labels': False, 'raw': False, 'extra_info': False, 'ignore_type': [], 'save_plot': False, 'save_loc': '../default/'}
        pltr.residual_plotter(path, ['z', 'Redshift'], **res_plt_args)
        pltr.residual_plotter(path, ['host_mass', 'Host Mass'], **res_plt_args)
    if 'hist' in types:
        hist_plt_arg = {'raw': False, 'save_plot': False, 'save_loc': '../default/', 'ignore_type': [], 'param_bins': [None, None, None, None, None, None]}
        pltr.snpy_histogram_plotter(path, **hist_plt_arg)
    return
def fit_SALT3(choice='CSP', single=''):
    salt_fit_params = {'mag_unc_max': 0.17, 'flux_unc_max': 625, 'processes': [],
                       'use_saved': False, 'show_plots': True, 'save': True, 'quiet': False}
    if len(single) > 0:
        salt_fit.indivisual_process(SN=single, data_set=choice, **snpy_fit_params)
    if choice in ['ATLAS', 'ZTF', 'CSP']:
        salt_fit.dataset_process(choice, **salt_fit_params)
    # elif choice in ['COMBINED']:
    #     salt_fit.combined_process(**snpy_fit_params)
    return
def plot_SALT3():
    path = CONSTANTS['salt_atlas_loc'] + 'atlas_salt_saved.txt'
    res_plt_args = {'sigma': [1, 1], 'labels': False, 'raw': False, 'extra_info': False, 'ignore_type': [],
                    'save_plot': False, 'save_loc': '../default/'}
    # pltr.residual_plotter(path, ['z', 'Redshift'], **plt_args)
    # pltr.residual_plotter(path, ['host_mass', 'Host Mass'], **plt_args)

    hist_plt_arg = {'raw': False, 'save_plot': False, 'save_loc': '../default/', 'ignore_type': [],
                    'param_bins': [None, None, None, None, None, None]}
    pltr.salt3_histogram_plotter(path, **hist_plt_arg)
    return

if __name__ == '__main__':
    start = systime.time() # Runtime tracker

    # fit_SNooPy('ATLAS')
    # plot_SNooPy('COMBINED')

    fit_SALT3('ATLAS')
    # plot_SALT3()

    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')








