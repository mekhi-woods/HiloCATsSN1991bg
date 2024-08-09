import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import numpy as np
    import snpy
    import glob
    import time as systime
    import scripts.general as gen
    import scripts.snpyfit as snpy_fit
    import scripts.saltfit as salt_fit
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
def plotter(choices, plot_type, cut=False, ignore_type=[]):
    for data_set in choices:
        path = CONSTANTS[data_set + '_saved_loc'] + data_set + '_saved.txt'
        if cut:
            path = CONSTANTS[data_set + '_saved_loc'] + data_set + '_saved_cut.txt'

        if 'resid' in plot_type:
            print("Plotting Hubble Residuals for '"+data_set+"'...")
            res_plot_args = {'path': path, 'save_loc': CONSTANTS[data_set+'_saved_loc'], 'sigma': [1, 1],
                             'labels': False, 'raw': False, 'extra_info': True, 'ignore_type': ignore_type, 'save_plot': True}
            gen.residual_plotter(x_params=['z', 'Redshift'], **res_plot_args)
            gen.residual_plotter(x_params=['host_mass', 'Host Mass'], **res_plot_args)

        if 'hist' in plot_type:
            print("Plotting parameter histograms for '"+data_set+"'...")
            hist_args = {'save_loc': CONSTANTS[data_set+'_saved_loc'], 'ignore_type': ignore_type,
                         'raw': False, 'save_plot': False, 'param_bins': [None, None, None, None, None]}
            gen.snpy_histogram_plotter(path=path, **hist_args)
    return

if __name__ == '__main__':
    start = systime.time()  # Runtime tracker

    # SNooPy Fitting
    snpy_fit_params = {'replot': False, 'use_saved': True, 'show_plots': True, 'save': False, 'quiet': False}
    snpy_fit.dataset_process(data_set='ATLAS', **snpy_fit_params) # 'CSP', 'ATLAS', 'ZTF'
    # snpy_fit.combined_process(**snpy_fit_params)
    # snpy_fit.indivisual_process(SN='2009F', source='ATLAS', **snpy_fit_params)

    # SALT3 Fitting
    salt_fit_params = {}
    # salt_fit.dataset_process('CSP', **salt_fit_params) # 'CSP', 'ATLAS', 'ZTF'

    # fit_reviewer('COMBINED')

    # ignore_type = []
    # gen.snpy_sample_cutter('csp_saved_loc')
    # *choices* : combined/atlas/ztf/csp || *plot_type* : resid/hist
    # plotter(['combined'], ['resid'], cut=True, ignore_type=ignore_type)

# -------------------------------------------------------------------------------------------------------------------- #

    # txt_file = '../saved/snpy/combined/combined_saved_cut.txt'
    # cuts, mass_step_arr, mass_step_err_arr = np.arange(8, 11, 0.1), np.array([]), np.array([])
    # for c in cuts:
    #     m, m_err = gen.mass_step_calc(txt_file, cut=c, quiet=True)
    #     mass_step_arr = np.append(mass_step_arr, m)
    #     mass_step_err_arr = np.append(mass_step_err_arr, m_err)
    #
    # max_step = np.max(mass_step_arr)
    # max_step_err = mass_step_err_arr[np.where(mass_step_arr == np.max(mass_step_arr))[0][0]]
    # gen.mass_step_calc(txt_file, cut=max_step, quiet=False)
    #
    # import matplotlib.pyplot as plt
    # plt.figure(figsize=(12, 6))
    # plt.title('Mass Step Cuts\n' + 'Max: ' + str(round(max_step, 2)) + ' +/- ' + str(round(max_step_err, 2)))
    # plt.errorbar(cuts, mass_step_arr, yerr=mass_step_err_arr, fmt='o')
    # plt.xlabel('Host Mass Cutting Point'); plt.ylabel('Mass Step')
    # plt.show()

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








