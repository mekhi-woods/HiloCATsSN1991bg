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

# def fit_reviewer(choice='COMBINED'):
#     if choice == 'COMBINED':
#         saved_models_loc = [glob.glob(CONSTANTS['csp_saved_loc'] + 'models/*.snpy'),
#                             glob.glob(CONSTANTS['combined_saved_loc'] + 'models/*.snpy')]
#     else:
#         raise ValueError("[!!!] No matches! ['COMBINED']")
#
#     all_bad_fit, all_okay_fit, all_good_fit = [], [], []
#     for models in saved_models_loc:
#         print('========================================================================')
#         print('Fits for', models[0].split('/')[-2].upper(), 'data...')
#         print('========================================================================')
#         bad_fit, okay_fit, good_fit = [], [], []
#         for model in models:
#             objname = model.split('/')[-1].split('_')[0]
#             if objname[:2] == 'SN':
#                 objname = objname[2:]
#             print('------------------------------------------------------------------------')
#             print('[', models.index(model)+1, '/', len(models), '] - ', objname)
#
#             pltr.lc_replot(model, spread=[30, 30], stacked=False)
#
#             selecting = True
#             while selecting:
#                 selection = input('How is this fit? [1-bad, 2-okay, 3-good]: ')
#                 if selection == '1':
#                     bad_fit.append(objname)
#                     selecting = False
#                 elif selection == '2':
#                     okay_fit.append(objname)
#                     selecting = False
#                 elif selection == '3':
#                     good_fit.append(objname)
#                     selecting = False
#                 else:
#                     print('[!!!] Unknown selection.')
#         print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
#         print('Good fits: ', good_fit)
#         print('Okay fits: ', okay_fit)
#         print('Bad fits: ', bad_fit)
#
#         all_good_fit = all_good_fit + good_fit
#         all_okay_fit = all_okay_fit + okay_fit
#         all_bad_fit = all_bad_fit + bad_fit
#
#
#     all_fits = {'good': all_good_fit, 'okay': all_okay_fit, 'bad': all_bad_fit}
#     with open(CONSTANTS['reviewed_fits_txt'], 'w') as f:
#         for fits in all_fits:
#             f.write(fits)
#             for obj in all_fits[fits]:
#                 f.write(', '+"'"+obj+"'")
#             f.write('\n')
#     return
# def plotter(choices, plot_type, cut=False, ignore_type=[]):
#     for data_set in choices:
#         path = CONSTANTS[data_set + '_saved_loc'] + data_set + '_saved.txt'
#         if cut:
#             path = CONSTANTS[data_set + '_saved_loc'] + data_set + '_saved_cut.txt'
#
#         if 'resid' in plot_type:
#             print("Plotting Hubble Residuals for '"+data_set+"'...")
#             res_plot_args = {'path': path, 'save_loc': CONSTANTS[data_set+'_saved_loc'], 'sigma': [1, 1],
#                              'labels': False, 'raw': False, 'extra_info': True, 'ignore_type': ignore_type, 'save_plot': True}
#             pltr.residual_plotter(x_params=['z', 'Redshift'], **res_plot_args)
#             pltr.residual_plotter(x_params=['host_mass', 'Host Mass'], **res_plot_args)
#
#         if 'hist' in plot_type:
#             print("Plotting parameter histograms for '"+data_set+"'...")
#             hist_args = {'save_loc': CONSTANTS[data_set+'_saved_loc'], 'ignore_type': ignore_type,
#                          'raw': False, 'save_plot': False, 'param_bins': [None, None, None, None, None]}
#             pltr.snpy_histogram_plotter(path=path, **hist_args)
#     return
def fit_SALT3():
    salt_fit_params = {'use_saved': False, 'replot': False, 'show_plots': True, 'save': True, 'quiet': False}
    salt_fit.dataset_process('ATLAS', **salt_fit_params) # 'CSP', 'ATLAS', 'ZTF'
    return
def fit_SNooPy():
    snpy_fit_params = {'mag_unc_max': 0, 'flux_unc_max': 0,
                       'replot': False, 'use_saved': True, 'show_plots': True, 'save': True, 'quiet': False}
    # snpy_fit.dataset_process(data_set='CSP', **snpy_fit_params) # 'CSP', 'ATLAS', 'ZTF'
    # snpy_fit.combined_process(**snpy_fit_params)
    snpy_fit.indivisual_process(SN='req00381099', data_set='ZTF', **snpy_fit_params)
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
    start = systime.time()  # Runtime tracker

    # fit_SNooPy()
    # fit_SALT3()
    # plot_SALT3()

    objs = gen.dict_handler(choice='arrays', path=CONSTANTS['salt_atlas_loc'] + 'atlas_salt_saved.txt')

    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')








