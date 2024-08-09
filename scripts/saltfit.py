import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

import scripts.general as gen

CONSTANTS = gen.get_constants()

def dataset_process(data_set):
    if data_set == 'CSP':
        # Set objects
        objs = gen.data_proccesser(data_set='CSP', quiet=True)

        # Creating tables for SALT3
        cspTbls = {}
        for obj in objs:
            cspTbls.update({obj: Table([objs[obj]['time'], objs[obj]['filters'], objs[obj]['flux'],
                                        objs[obj]['dflux'], objs[obj]['zp'],
                                        np.full(len(objs[obj]['time']), 'ab')],
                                        names=('time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys')) })




        # filter_wheel = {'B': 'B', 'H': 'H', 'J': 'J', 'Jrc2': 'Jrc2', 'V': 'V', 'Y': 'Y', 'Ydw': 'Ydw', 'g': 'g',
        #                 'i': 'i', 'r': 'r', 'u': 'u'}
        # gen.write_ASCII(objs=cspObjs, filter_set=filter_wheel, save_loc=CONSTANTS['csp_saved_loc']+'ascii/')
        #
        # # Replot
        # if replot:
        #     plt_args = {'color_wheel': [None, None, None, None, None, None, None, None, None, None, None],
        #                 'save_loc': CONSTANTS['csp_saved_loc']+'plots/', 'y_type': 'mag', 'pause_time': 2,
        #                 'quiet': quiet, 'save_plots': save}
        #     gen.lc_plot(objs=cspObjs, **plt_args)
        #
        # # Fit with SNooPy
        # fit_args = {'save_loc': CONSTANTS['csp_saved_loc'], 'use_saved': False, 'snpy_plots': True, 'save_plots': save, 'quiet': quiet}
        # cspParams = gen.snpy_fit(paths=glob.glob(CONSTANTS['csp_saved_loc']+'ascii/*.txt'), **fit_args)
        #
        # # Save parameters to file
        # gen.dict_handler(choice='pack', data_dict=cspParams, path=CONSTANTS['csp_saved_loc']+'csp_saved.txt')
        #
        # # Get host masses
        # gen.host_mass(CONSTANTS['csp_saved_loc']+'csp_saved.txt', save_loc='../default/', keep_data=False, update_saved=True)
    elif data_set == 'ATLAS':
        x=1
    elif data_set == 'ZTF':
        x=1
    else:
        raise ValueError("Data set not supported ['CSP'/'ATLAS'/'ZTF']")
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