# import os
# import sys
# import snpy
# import shutil
# import time as systime

# from zipfile import ZipFile
# from scripts import tns_redshifts
# from astropy import cosmology as cosmo
# from astropy.table import QTable, Table, Column
#
# from astro_ghost.ghostHelperFunctions import getTransientHosts
# from astropy.coordinates import SkyCoord
# from astropy.cosmology import FlatLambdaCDM
# from astropy import units as u
# import sncosmo
#
# from astroquery.sdss import SDSS
# from astroquery.mast import Catalogs

import glob
import numpy as np
import matplotlib.pyplot as plt

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
    gen.write_ASCII(objs=objs, **ASCIIArgs) # Write ASCII files for SNooPy fitting
    if replot:
        gen.lc_plot(objs=objs, **plottingArgs) # Replot LCs
    objParams = gen.snpy_fit(paths=glob.glob(save_loc + 'ascii/*.txt'), **fittingArgs) # Fit with SNooPy
    gen.dict_handler(choice='pack', data_dict=objParams, path=save_txt) # Save parameters to file
    gen.host_mass(save_txt, keep_data=False, update_saved=True) # Get host masses

    return
# def combined_process(replot=True, save=False, quiet=False):
#     save_loc = CONSTANTS['combined_saved_loc']
#     save_txt = save_loc + 'combined_saved.txt'
#     processingArgs = {'data_set': data_set, 'quiet': True}
#     fittingArgs = {'save_loc': save_loc, 'use_saved': False, 'snpy_plots': True, 'save_plots': save, 'quiet': quiet}
#     plottingArgs = {'save_loc': save_loc + 'plots/', 'y_type': 'mag', 'pause_time': 2, 'quiet': quiet, 'save_plots': save}
#     ASCIIArgs = {'filter_set': {'c': 'ATgr', 'o': 'ATri', 'ZTF_g': 'g', 'ZTF_r': 'r', 'ZTF_i': 'i'}, 'save_loc': save_loc + 'ascii/'}
#
#     objs = gen.alt_data_proccesser(**processingArgs) # Set objects
#     gen.write_ASCII(objs=objs, **ASCIIArgs) # Write ASCII files for SNooPy fitting
#     if replot:
#         gen.lc_plot(objs=objs, **plottingArgs) # Replot LCs
#     objParams = gen.snpy_fit(paths=glob.glob(save_loc + 'ascii/*.txt'), **fittingArgs) # Fit with SNooPy
#     gen.dict_handler(choice='pack', data_dict=objParams, path=save_txt) # Save parameters to file
#     gen.host_mass(save_txt, keep_data=False, update_saved=True) # Get host masses
#     return
# def indivisual_process(SN, source, replot=True, save=False, quiet=False):
#
#     if source == 'CSP':
#         f_w = {'B': 'B', 'H': 'H', 'J': 'J', 'Jrc2': 'Jrc2', 'V': 'V', 'Y': 'Y',
#                'Ydw': 'Ydw', 'g': 'g', 'i': 'i', 'r': 'r', 'u': 'u'}
#         color_wheel = {'B': None, 'H': None, 'J': None, 'Jrc2': None, 'V': None, 'Y': None,
#                        'Ydw': None, 'g': None, 'i': None, 'r': None, 'u': None}
#         ascii_args = {'filter_set': f_w, 'save_loc': CONSTANTS['csp_saved_loc']+'ascii/'}
#         plt_args = {'color_wheel': [None, None, None, None, None, None, None, None, None, None, None],
#                     'save_loc': CONSTANTS['csp_saved_loc'] + 'plots/', 'y_type': 'mag', 'pause_time': 2,
#                     'quiet': quiet, 'save_plots': save}
#         fit_args = {'save_loc': CONSTANTS['csp_saved_loc'], 'use_saved': False, 'snpy_plots': True,
#                     'save_plots': save, 'quiet': quiet}
#         save_args = {'path': CONSTANTS['csp_saved_loc']+'csp_saved_indv.txt'}
#         replot_args = {'save_plot': False, 'stacked': False, 'save_loc': '../default/',
#                        'colors': {'g': 'green', 'r': 'red', 'i': 'violet'}, 'spread': [30, 30]}
#     elif source == 'ATLAS':
#         ascii_args = {'filter_set': {'c': 'ATgr', 'o': 'ATri'},
#                       'save_loc': CONSTANTS['atlas_saved_loc']+'ascii/'}
#         plt_args = {'color_wheel': ['orange', 'cyan', 'violet'],
#                     'save_loc': CONSTANTS['atlas_saved_loc'] + 'plots/', 'y_type': 'mag', 'pause_time': 2,
#                     'quiet': False, 'save_plots': True}
#         fit_args = {'save_loc': CONSTANTS['atlas_saved_loc'], 'use_saved': False, 'snpy_plots': True,
#                     'save_plots': True, 'quiet': False}
#         save_args = {'path': CONSTANTS['ztf_saved_loc']+'ztf_saved_indv.txt'}
#         replot_args = {'save_plot': False, 'stacked': False, 'save_loc': '../default/',
#                        'colors': {'g': 'green', 'r': 'red', 'i': 'violet'}, 'spread': [30, 30]}
#     elif source == 'ZTF':
#         ascii_args = {'filter_set': {'ZTF_g': 'g', 'ZTF_r': 'r', 'ZTF_i': 'i'},
#                       'save_loc': CONSTANTS['ztf_saved_loc']+'ascii/'}
#         plt_args = {'color_wheel': ['green', 'red', 'violet', 'blue'],
#                     'save_loc': CONSTANTS['ztf_saved_loc'] + 'plots/', 'y_type': 'mag', 'pause_time': 1,
#                     'quiet': False, 'save_plots': True}
#         fit_args = {'paths': ['../saved/ztf/ascii/'+target_SN+'_snpy.txt'], 'save_loc': CONSTANTS['atlas_saved_loc'],
#                     'use_saved': True, 'snpy_plots': True, 'save_plots': True, 'quiet': False}
#         save_args = {'path': CONSTANTS['ztf_saved_loc']+'ztf_saved_indv.txt'}
#         replot_args = {'save_plot': False, 'stacked': False, 'save_loc': '../default/',
#                        'colors': {'g': 'green', 'r': 'red', 'i': 'violet'}, 'spread': [30, 30]}
#     else:
#         raise ValueError("Data set must be either 'ATLAS' or 'ZTF'")
#     objs = gen.alt_data_proccesser(data_set=data_set, quiet=True)
#     gen.write_ASCII(objs={target_SN: objs[target_SN]}, **ascii_args)
#     gen.lc_plot(objs={target_SN: objs[target_SN]} , **plt_args)
#     params = gen.snpy_fit(**fit_args)
#     gen.dict_handler(choice='pack', data_dict=params, **save_args)
#     # gen.lc_replot(CONSTANTS['ztf_saved_loc']+'models/'+target_SN+'_EBV_model2.snpy', **replot_args)
#     return
# def alt_indivisual_process(SN, source, replot=True, save=False, quiet=False):
#     save_loc = CONSTANTS[data_set.lower() + '_saved_loc']
#     save_txt = save_loc + data_set.lower() + '_saved.txt'
#     processingArgs = {'data_set': data_set, 'quiet': True}
#     fittingArgs = {'save_loc': save_loc, 'use_saved': False, 'snpy_plots': True, 'save_plots': save, 'quiet': quiet}
#     plottingArgs = {'save_loc': save_loc + 'plots/', 'y_type': 'mag', 'pause_time': 2, 'quiet': quiet, 'save_plots': save}
#
#     if data_set == 'CSP':
#         ASCIIArgs = {'filter_set': {'B': 'B', 'H': 'H', 'J': 'J', 'Jrc2': 'Jrc2', 'V': 'V', 'Y': 'Y', 'Ydw': 'Ydw',
#                                     'g': 'g', 'i': 'i', 'r': 'r', 'u': 'u'}, 'save_loc': save_loc + 'ascii/'}
#     elif data_set == 'ATLAS':
#         ASCIIArgs = {'filter_set': {'c': 'ATgr', 'o': 'ATri'},
#                      'save_loc': save_loc + 'ascii/'}
#     elif data_set == 'ZTF':
#         ASCIIArgs = {'filter_set': {'ZTF_g': 'g', 'ZTF_r': 'r', 'ZTF_i': 'i'},
#                      'save_loc': save_loc + 'ascii/'}
#     else:
#         raise ValueError("Data set not supported ['CSP'/'ATLAS'/'ZTF']")
#
#     objs = gen.alt_data_proccesser(**processingArgs) # Set objects
#     gen.write_ASCII(objs=objs, **ASCIIArgs) # Write ASCII files for SNooPy fitting
#     if replot:
#         gen.lc_plot(objs=objs, **plottingArgs) # Replot LCs
#     objParams = gen.snpy_fit(paths=glob.glob(save_loc + 'ascii/*.txt'), **fittingArgs) # Fit with SNooPy
#     gen.dict_handler(choice='pack', data_dict=objParams, path=save_txt) # Save parameters to file
#     gen.host_mass(save_txt, keep_data=False, update_saved=True) # Get host masses
#
#     return
# def write_ASCII(objs, filter_set, save_loc):
#     print('[+++] Saving data to ASCII files for SNooPy...')
#     for obj in objs:
#         with open(save_loc + obj + '_snpy.txt', 'w') as f:
#             # Line 1 -- Objname, Helio-Z, RA, Dec (Ex. SN1981D 0.005871 50.65992 -37.23272)
#             f.write(str(obj)+' '+str(objs[obj]['z'])+' '+str(objs[obj]['ra'])+' '+str(objs[obj]['dec'])+'\n')
#             for f_w in filter_set:
#                 f_indexs = np.where(objs[obj]['filters'] == f_w)[0]
#                 f.write('filter '+filter_set[f_w]+'\n')
#                 for i in f_indexs:
#                     # filter photometry block -- Date (JD/MJD), mag, err (i.e. 674.8593 12.94 0.11)
#                     f.write(str(objs[obj]['time'][i])+'\t'+str(objs[obj]['mag'][i])+'\t'+str(objs[obj]['dmag'][i])+'\n')
#     return
# def snpy_fit(paths, save_loc, use_saved=False, snpy_plots=True, save_plots=True, quiet=False):
#     print('[+++] Fitting data with SNooPy...')
#
#     # Check quiet
#     if quiet:
#         sys.stdout = open(os.devnull, 'w')
#
#     objParams = {}
#     for path in paths:
#         objname = path.split('/')[-1][:-9]
#         print('-------------------------------------------------------------------------------------------------------')
#         print('[', paths.index(path)+1, '/', len(paths), '] --', objname)
#
#         # Check saved fits
#         saved_found = False
#         if use_saved:
#             for saved_fit in glob.glob(save_loc+'/models/*.snpy'):
#                 if objname == saved_fit.split('/')[-1][:-16]:
#                     print('Saved fit found for '+objname+'...')
#                     saved_found = True
#                     n_s = snpy.get_sn(saved_fit)
#                     mjds, mjde = 999999999, 0
#                     for filter in list(n_s.data.keys()):
#                         n_min, n_max = min(n_s.data[filter].MJD), max(n_s.data[filter].MJD)
#                         if mjds > n_min:
#                             mjds = n_min
#                         if mjde < n_max:
#                             mjde = n_max
#                     objParams.update({objname: {'ra': n_s.ra, 'dec': n_s.decl, 'z': n_s.z, 'MJDs': mjds, 'MJDe': mjde,
#                                                 'mu': n_s.parameters['DM'], 'st': n_s.parameters['st'],
#                                                 'Tmax': n_s.parameters['Tmax'], 'EBVhost': n_s.parameters['EBVhost'],
#                                                 'mu_err': n_s.errors['DM'], 'st_err': n_s.errors['st'],
#                                                 'Tmax_err': n_s.errors['Tmax'], 'EBVhost_err': n_s.errors['EBVhost']}})
#                     if snpy_plots:
#                         n_s.plot()
#                         plt.show()
#                         systime.sleep(3)
#                     print(objParams[objname])
#         if saved_found:
#             continue
#
#         problem_children = handle_problem_children(state='READ') # Open problem children
#
#         # Load data
#         try:
#             n_s = snpy.get_sn(path)
#             n_s.choose_model('EBV_model2', stype='st')
#             n_s.set_restbands()  # Auto pick appropriate rest-bands
#         except:
#             print('[!!!] Failed to load ASCII file')
#             continue
#
#         # Remove empty filters -- fix for 'ValueError: attempt to get argmin of an empty sequence'
#         for filter in list(n_s.data.keys()):
#             if len(n_s.data[filter].magnitude) == 0:
#                 del n_s.data[filter]
#         print('Best filters:', list(n_s.data.keys()))
#
#         run = True
#         while run:
#             try:
#                 # Minimum MJD and Max MJD
#                 mjds, mjde = 999999999, 0
#                 for filter in list(n_s.data.keys()):
#                     n_min, n_max = min(n_s.data[filter].MJD), max(n_s.data[filter].MJD)
#                     if mjds > n_min:
#                         mjds = n_min
#                     if mjde < n_max:
#                         mjde = n_max
#
#                 n_s.k_version = '91bg'
#                 n_s.fit(bands=None, dokcorr=True, k_stretch=False, reset_kcorrs=True, **{'mangle': 1, 'calibration': 0})
#                 n_s.save(save_loc + 'models/' + objname + '_EBV_model2.snpy')
#                 run = False
#
#                 if snpy_plots:
#                     n_s.plot(outfile=save_loc + 'plots/' + objname + '_snpyplots.png')
#                     plt.show()
#                 plt.close()
#
#                 objParams.update({objname : {'ra': n_s.ra, 'dec': n_s.decl, 'z': n_s.z, 'MJDs': mjds, 'MJDe': mjde,
#                                              'mu': n_s.parameters['DM'], 'st': n_s.parameters['st'],
#                                              'Tmax': n_s.parameters['Tmax'], 'EBVhost': n_s.parameters['EBVhost'],
#                                              'mu_err': n_s.errors['DM'], 'st_err': n_s.errors['st'],
#                                              'Tmax_err': n_s.errors['Tmax'], 'EBVhost_err': n_s.errors['EBVhost']}})
#             except Exception as error:
#                 if 'All weights for filter' and 'are zero.' in str(error):
#                     print('[!!!] Weights for filter', str(error).split(' ')[4], 'are zero. Removing...')
#                     del n_s.data[str(error).split(' ')[4]]
#                 elif 'Error:  to solve for EBVhost, you need to fit more than one filter' in str(error):
#                     print('[!!!]', error)
#                     problem_children = np.append(problem_children, objname)
#                     handle_problem_children(state='WRITE', problem_c=problem_children) # Commit problem children
#                     run = False
#                 else:
#                     print(error)
#                     problem_children = np.append(problem_children, objname)
#                     handle_problem_children(state='WRITE', problem_c=problem_children) # Commit problem children
#                     run = False
#
#     print('Successfully fit [', len(objParams), '/', len(paths), '] !')
#
#     # Restore print statements
#     sys.stdout = sys.__stdout__
#
#     return objParams
# def snpy_sample_cutter(save_loc):
#     print('[+++] Cutting sample for ' + save_loc.split('_')[0].upper() + ' data...')
#
#     path = glob.glob(get_constants()[save_loc] + '*_saved.txt')[0]
#     objs = dict_handler(choice='unpack', path=path)
#     i = 0
#     new_objs = {}
#     for obj in objs:
#         print('[' + str(list(objs.keys()).index(obj) + 1) + '/' + str(len(objs)) + '] -- ' + obj)
#         print('---------------------------------------------------------------------------------------------------')
#         resid = float(objs[obj]['mu']) - CURRENT_COSMO.distmod(float(objs[obj]['z'])).value
#         resid -= np.median(resid)
#
#         if float(objs[obj]['EBVhost']) < -0.2 or float(objs[obj]['EBVhost']) > 0.2:
#             print('[!!!] EBVhost out of range.')
#             continue
#         if float(objs[obj]['EBVhost_err']) > 0.1:
#             print('[!!!] EBVhost errors out of range.')
#             continue
#
#         if float(objs[obj]['st']) < 0.3 or float(objs[obj]['st']) > 1.0:
#             print('[!!!] Stretch out of range.')
#             continue
#         if float(objs[obj]['st_err']) > 0.1:
#             print('[!!!] Stretch error out of range.')
#             continue
#
#         if float(objs[obj]['Tmax_err']) > 1:
#             print('[!!!] Maximum time error out of range.')
#             continue
#
#         if float(objs[obj]['z']) < 0.015:
#             print('[!!!] Redshift out of range.')
#             continue
#
#         if float(objs[obj]['host_mass']) == 0.00:
#             print('[!!!] Host mass null.')
#             continue
#         i = i + 1
#
#         # Save obj to new dict
#         new_objs.update({obj: objs[obj]})
#
#     dict_handler(data_dict=new_objs, choice='pack', path=get_constants()[save_loc] + save_loc.split('_')[0] + '_saved_cut.txt')
#     return
