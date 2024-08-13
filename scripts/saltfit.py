import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

import scripts.general as gen

CONSTANTS = gen.get_constants()

def dataset_process(data_set, use_saved=True, replot=True, show_plots=True, save=False, quiet=False):
    save_loc = CONSTANTS[data_set.lower() + '_saved_loc']
    save_txt = save_loc + data_set.lower() + '_saved.txt'
    processingArgs = {'data_set': data_set, 'quiet': True}
    fittingArgs = {'save_loc': save_loc, 'use_saved': use_saved, 'snpy_plots': show_plots, 'save_plots': save, 'quiet': quiet}
    plottingArgs = {'save_loc': save_loc + 'plots/', 'y_type': 'mag', 'pause_time': 2, 'quiet': quiet, 'save_plots': save}

    if data_set == 'CSP':
        ASCIIArgs = {'filter_set': {'B': 'B', 'H': 'H', 'J': 'J', 'Jrc2': 'Jrc2', 'V': 'V', 'Y': 'Y', 'Ydw': 'Ydw',
                                    'g': 'g', 'i': 'i', 'r': 'r', 'u': 'u'}, 'save_loc': save_loc + 'ascii/'}
    elif data_set == 'ATLAS':
        ASCIIArgs = {'filter_set': {'c': 'ATgr', 'o': 'ATri'},
                     'save_loc': save_loc + 'ascii/'}
    elif data_set == 'ZTF':
        ASCIIArgs = {'filter_set': {'ZTF_g': 'g', 'ZTF_r': 'r', 'ZTF_i': 'i'},
                     'save_loc': save_loc + 'ascii/'}
    else:
        raise ValueError("Data set not supported ['CSP'/'ATLAS'/'ZTF']")

    objs = gen.data_proccesser(**processingArgs) # Set objects
    write_ASCII(objs=objs, **ASCIIArgs) # Write ASCII files for SNooPy fitting
    if replot:
        gen.lc_plot(objs=objs, **plottingArgs) # Replot LCs
    objParams = snpy_fit(paths=glob.glob(save_loc + 'ascii/*.txt'), **fittingArgs) # Fit with SNooPy
    gen.dict_handler(choice='pack', data_dict=objParams, path=save_txt) # Save parameters to file
    gen.host_mass(save_txt, keep_data=False, update_saved=True) # Get host masses

    return
def salt3_atlas_process():
    objs = gen.data_proccesser(data_set='ATLAS', mag_unc_max=0, flux_unc_max=0, quiet=True)
    for obj in objs:
        if len(objs[obj]['filters']) == 0:
            continue
        objs[obj]['filters'][np.where(objs[obj]['filters'] == 'c')[0]] = 'atlasc'
        objs[obj]['filters'][np.where(objs[obj]['filters'] == 'o')[0]] = 'atlaso'

    # gen.lc_plot(objs, y_type='flux', pause_time=1, color_wheel=['orange', 'cyan', 'violet'], quiet=False, save_plots=False)

        # break


    # fit_args = {'plot_data': True, 'save_plot': False}
    # atlasParams = gen.salt3_fit(objs, plot_save_loc=gen.get_constants()['salt_atlas_loc'], **fit_args)
    # gen.dict_handler(data_dict=atlasParams, choice='pack', path=CONSTANTS['salt_atlas_loc']+'salt_atlas_saved.txt')

    # gen.host_mass(CONSTANTS['salt_atlas_loc']+'salt_atlas_saved.txt', save_loc='../default/', keep_data=False, update_saved=True)
    #
    # res_plot_args = {'save_loc': CONSTANTS['salt_atlas_loc'], 'sigma': [1, 1],
    #                  'labels': False, 'raw': False, 'extra_info': True, 'ignore_type': [], 'save_plot': True}
    # gen.residual_plotter(path=CONSTANTS['salt_atlas_loc']+'salt_atlas_saved.txt', x_params=['z', 'Redshift'], **res_plot_args)
    # gen.residual_plotter(path=CONSTANTS['salt_atlas_loc']+'salt_atlas_saved.txt', x_params=['host_mass', 'Host Mass'], **res_plot_args)

    return
def salt3_ztf_process():
    objs = gen.data_proccesser(data_set='ZTF', mag_unc_max=0, flux_unc_max=0, quiet=False)
    # Swap filters to SALT3 names
    for obj in objs:
        if len(objs[obj]['filters']) == 0:
            continue
        objs[obj]['filters'][np.where(objs[obj]['filters'] == 'ZTF_g')[0]] = 'ztfg'
        objs[obj]['filters'][np.where(objs[obj]['filters'] == 'ZTF_r')[0]] = 'ztfr'
        objs[obj]['filters'][np.where(objs[obj]['filters'] == 'ZTF_i')[0]] = 'ztfi'

    # gen.lc_plot(objs, y_type='flux', pause_time=1, color_wheel=['green', 'red', 'violet'], quiet=False, save_plots=False)

    fit_args = {'plot_save_loc': CONSTANTS['salt_ztf_loc']+'plots/', 'plot_data': True, 'save_plot': True}
    atlasParams = gen.salt3_fit(objs, **fit_args)
    gen.dict_handler(data_dict=atlasParams, choice='pack', path=CONSTANTS['salt_ztf_loc']+'salt_ztf_saved.txt')

    # gen.host_mass(CONSTANTS['salt_ztf_loc']+'salt_ztf_saved.txt', save_loc='../default/', keep_data=False, update_saved=True)
    #
    # res_plot_args = {'save_loc': CONSTANTS['salt_ztf_loc'], 'sigma': [1, 1],
    #                  'labels': False, 'raw': False, 'extra_info': True, 'ignore_type': [], 'save_plot': True}
    # gen.residual_plotter(path=CONSTANTS['salt_ztf_loc']+'salt_ztf_saved.txt', x_params=['z', 'Redshift'], **res_plot_args)
    # gen.residual_plotter(path=CONSTANTS['salt_ztf_loc']+'salt_ztf_saved.txt', x_params=['host_mass', 'Host Mass'], **res_plot_args)

    return