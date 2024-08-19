import glob
import sncosmo
import numpy as np
import time as systime
import matplotlib.pyplot as plt
from astropy.table import Table

import scripts.general as gen
import scripts.plotter as pltr

CONSTANTS = gen.get_constants()

def dataset_process(data_set, use_saved=True, replot=True, show_plots=True, save=False, quiet=False):
    save_loc = CONSTANTS['salt_'+data_set.lower()+'_loc']
    save_txt = save_loc + data_set.lower() + '_salt_saved.txt'
    processingArgs = {'data_set': data_set, 'quiet': False}
    fittingArgs = {'plot_save_loc': '../default/', 'plot_data': True, 'save_plot': True}
    plottingArgs = {'save_loc': save_loc + 'plots/', 'y_type': 'mag', 'pause_time': 2, 'quiet': quiet, 'save_plots': save}

    if data_set not in ['CSP', 'ATLAS', 'ZTF']:
        raise ValueError("Data set not supported ['CSP'/'ATLAS'/'ZTF']")

    objs = gen.data_proccesser(**processingArgs) # Set objects

    # Fix filters for SALT3
    for obj in objs:
        if len(objs[obj]['filters']) == 0:
            continue
        if data_set == 'CSP':
            objs[obj]['filters'][np.where(objs[obj]['filters'] == 'u')[0]] = 'cspu'
            objs[obj]['filters'][np.where(objs[obj]['filters'] == 'g')[0]] = 'cspg'
            objs[obj]['filters'][np.where(objs[obj]['filters'] == 'r')[0]] = 'cspr'
            objs[obj]['filters'][np.where(objs[obj]['filters'] == 'i')[0]] = 'cspi'
            objs[obj]['filters'][np.where(objs[obj]['filters'] == 'B')[0]] = 'cspB'
            objs[obj]['filters'][np.where(objs[obj]['filters'] == 'V0')[0]] = 'cspv3014'
            objs[obj]['filters'][np.where(objs[obj]['filters'] == 'V1')[0]] = 'cspv3009'
            objs[obj]['filters'][np.where(objs[obj]['filters'] == 'V')[0]] = 'cspv9844'
            objs[obj]['filters'][np.where(objs[obj]['filters'] == 'Y')[0]] = 'cspys'
            objs[obj]['filters'][np.where(objs[obj]['filters'] == 'J')[0]] = 'cspjs'
            objs[obj]['filters'][np.where(objs[obj]['filters'] == 'Jrc2')[0]] = 'cspjrc2'
            objs[obj]['filters'][np.where(objs[obj]['filters'] == 'Jdw')[0]] = 'cspjd'
            objs[obj]['filters'][np.where(objs[obj]['filters'] == 'Ydw')[0]] = 'cspyd'
            objs[obj]['filters'][np.where(objs[obj]['filters'] == 'Hdw')[0]] = 'csphd'
            objs[obj]['filters'][np.where(objs[obj]['filters'] == 'H')[0]] = 'csphs'
        elif data_set == 'ATLAS':
            objs[obj]['filters'][np.where(objs[obj]['filters'] == 'c')[0]] = 'atlasc'
            objs[obj]['filters'][np.where(objs[obj]['filters'] == 'o')[0]] = 'atlaso'
        elif data_set == 'ZTF':
            objs[obj]['filters'][np.where(objs[obj]['filters'] == 'ZTF_g')[0]] = 'ztfg'
            objs[obj]['filters'][np.where(objs[obj]['filters'] == 'ZTF_r')[0]] = 'ztfr'
            objs[obj]['filters'][np.where(objs[obj]['filters'] == 'ZTF_i')[0]] = 'ztfi'
        else:
            raise ValueError("Data set not supported ['CSP'/'ATLAS'/'ZTF']")
    if replot:
        pltr.lc_plot(objs=objs, **plottingArgs) # Replot LCs
    objParams = salt3_fit(objs=objs, **fittingArgs) # Fit with SALT3
    # objParams = salt3_sample_cutter(objParams) # Cutting
    gen.dict_handler(data_dict=objParams, choice='pack', path=save_txt) # Save data to file
    gen.host_mass(save_txt, save_loc='../default/', keep_data=False, update_saved=True) # Add masses to file
    return
def salt3_fit(objs, plot_save_loc='../default/', plot_data=True, save_plot=True):
    print('[+++] Fitting data with SALT3...')

    alpha, beta = float(CONSTANTS['salt_alpha']), float(CONSTANTS['salt_beta'])
    mB_const, M0 = float(CONSTANTS['salt_mB_const']), float(CONSTANTS['salt_absolute_mag'])

    params = {}
    for obj in objs:
        print('-------------------------------------------------------------------------------------------------------')
        print('[', list(objs.keys()).index(obj)+1, '/', len(objs), ']')

        try:
            data = Table([objs[obj]['time'], objs[obj]['filters'], objs[obj]['flux'], objs[obj]['dflux'],
                         objs[obj]['zp'], np.full(len(objs[obj]['time']), 'ab')],
                         names=('time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'))

            # Create Model
            model = sncosmo.Model(source='salt3')

            # Fit data to model
            model.set(z=objs[obj]['z'])  # set the model's redshift.
            result, fitted_model = sncosmo.fit_lc(data, model, ['t0', 'x0', 'x1', 'c'], bounds={'x1': (-5, 5)})

            # Save parameters
            params.update({obj: {'ra': objs[obj]['ra'], 'dec': objs[obj]['dec'], 'z': objs[obj]['z']}})
            for i in range(len(result.vparam_names)):
                n_param = result.param_names[i+1]
                params[obj].update({n_param: result.parameters[i+1]})
                params[obj].update({n_param+'_err': result.errors[n_param]})

            # Calculate
            pho_mB = -2.5*np.log10(params[obj]['x0']) + mB_const
            pho_mB_err = np.abs((-2.5*params[obj]['x0_err']) / (params[obj]['x0_err'] * np.log(10)))

            mu = pho_mB + (alpha*params[obj]['x1']) - (beta*params[obj]['c']) - M0
            mu_err = np.sqrt(pho_mB_err**2 + (np.abs(alpha)*params[obj]['x1_err'])**2 + (np.abs(alpha)*params[obj]['c_err'])**2)

            params[obj].update({'mu': mu, 'mu_err': mu_err})

            # Print Results
            print(obj, '|', '('+str(params[obj]['ra'])+', '+str(params[obj]['dec'])+'), z =', params[obj]['z'])
            for p in result.vparam_names:
                print(p, '|', params[obj][p], '+/-', params[obj][p+'_err'])
            print('mu', '|', params[obj]['mu'], '+/-', params[obj]['mu_err'])

            # Plot data with fit
            if plot_data:
                sncosmo.plot_lc(data, model=fitted_model, errors=result.errors, yfigsize=8)
                if save_plot:
                    print('Saving plots to', plot_save_loc+obj+'_salt3lc.png')
                    plt.savefig(plot_save_loc+obj+'_salt3lc.png')
                plt.show()

            print('Pausing for 1 seconds...')
            systime.sleep(1)
        except Exception as error:
            print(error)

    print('Successfully fit [', len(params), '/', len(objs), '] !')

    return params
def salt3_sample_cutter(objs):
    print('[+++] Cutting SALT3 results...')

    new_objs = {}
    for obj in objs:
        print('[' + str(list(objs.keys()).index(obj) + 1) + '/' + str(len(objs)) + '] -- ' + obj)
        print('---------------------------------------------------------------------------------------------------')

        if float(objs[obj]['x1']) < -3 or float(objs[obj]['x1']) > 3:
            print('[!!!] X1 out of range.')
            continue
        if float(objs[obj]['x1_err']) > 1:
            print('[!!!] X1 error too large.')
            continue

        if float(objs[obj]['c']) < -0.3 or float(objs[obj]['c']) > 0.3:
            print('[!!!] c out of range.')
            continue
        if float(objs[obj]['c_err']) > 0.1:
            print('[!!!] C error too large.')
            continue

        if float(objs[obj]['t0_err']) > 2:
            print('[!!!] T0 error too large.')
            continue

    # Save obj to new dict
        new_objs.update({obj: objs[obj]})

    print('Successfully cut data [', len(new_objs), '/', len(objs), '] !')

    return new_objs