import os
import sys
import snpy
import glob
import shutil
import time as systime
import numpy as np
import matplotlib.pyplot as plt
from zipfile import ZipFile
from scripts import tns_redshifts
from astropy import cosmology as cosmo
from astropy.table import QTable, Table, Column

from astro_ghost.ghostHelperFunctions import getTransientHosts
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import sncosmo

from astroquery.sdss import SDSS
from astroquery.mast import Catalogs

CURRENT_COSMO = FlatLambdaCDM(70, 0.3)  # Hubble Constant, Omega-Matter

def get_constants(cosntant_loc='../txts/constants.txt'):
    CONSTANTS = {}
    with open(cosntant_loc, 'r') as f:
        temp = f.readlines()
        for line in temp:
            line = line[:-1].split(', ')
            if len(line) != 2:
                continue
            CONSTANTS.update({line[0]: line[1]})
    return CONSTANTS
def TNS_details(ra, dec):
    tns_bot_id, tns_bot_name, tns_bot_api_key = '73181', 'YSE_Bot1', '0d771345fa6b876a5bb99cd5042ab8b5ae91fc67'

    # Code from David
    headers = tns_redshifts.build_tns_header(tns_bot_id, tns_bot_name)
    tns_api_url = f"https://www.wis-tns.org/api/get"

    # get the API URLs
    search_tns_url = tns_redshifts.build_tns_url(tns_api_url, mode="search")
    get_tns_url = tns_redshifts.build_tns_url(tns_api_url, mode="get")

    search_data = tns_redshifts.build_tns_search_query_data(tns_bot_api_key, ra, dec)
    transients = tns_redshifts.rate_limit_query_tns(search_data, headers, search_tns_url)

    get_data = tns_redshifts.build_tns_get_query_data(tns_bot_api_key, transients[0])
    transient_detail = tns_redshifts.rate_limit_query_tns(get_data, headers, get_tns_url)

    return transient_detail
def TNS_objname_z(obj_ra, obj_dec):
    print('Searching TNS for ['+str(obj_ra)+', '+str(obj_dec)+']...')
    tns_key = np.genfromtxt(get_constants()['tns_key_txt'], dtype=str, delimiter=', ', skip_header=1)

    run, num_tries = True, 10
    while run:
        try:
            obj_name, obj_z = '', 0.00

            # Check TNS key
            if len(np.shape(tns_key)) == 2:
                print('Checking TNS key...')
                objnames, ra, dec, z = tns_key[:, 2], tns_key[:, 0].astype(float), tns_key[:, 1].astype(float), tns_key[:, 3]
                for n in range(len(objnames)):
                    if abs(obj_ra - ra[n]) < 0.01:
                        obj_name, obj_z = objnames[n], z[n]
                        return obj_name, obj_z

            # Query TNS key
            print('Querying TNS...')
            details = TNS_details(obj_ra, obj_dec)
            obj_name, obj_z = details['objname'], details['redshift']
            with open(get_constants()['tns_key_txt'], 'a') as f:
                f.write(str(obj_ra) + ', ' + str(obj_dec) + ', ' + obj_name + ', ' + str(obj_z) + '\n') # Save to key
            return obj_name, obj_z
        except Exception as error:
            print('***********************************************************************************************')
            print(error)
            print('TNS timed out, pausing for 5 seconds...')
            print(num_tries, 'tries left...')
            print('***********************************************************************************************')
            systime.sleep(5)
            if num_tries == 0:
                raise RuntimeError(f'TNS completly timed out after {num_tries} tries.')
            else:
                num_tries -= 1
def dict_handler(data_dict={}, choice='', path='../default/dict.txt', delimiter=', '):
    if choice == 'unpack':
        print('[+++] Unpacking objects from '+path+'...')
        with open(path, 'r') as f:
            hdr = f.readline()[:-1].split(delimiter)
        data = np.genfromtxt(path, delimiter=delimiter, dtype=str, skip_header=1)
        if len(data) == 0:
            return {}
        temp_objs = {}
        for i in range(len(data[:, 0])):
            obj = data[:, 0][i]
            temp_objs.update({obj: {}})
            for j in range(len(hdr)):
                temp_objs[obj].update({hdr[j]: data[i, j]})
        return temp_objs
    elif choice == 'pack':
        print('[+++] Packing objects into text file...')
        if len(data_dict) == 0:
            print('[!!!] Attempted to pack an empty dictionary!')
            return
        catagories = list(data_dict[list(data_dict.keys())[0]].keys())
        with open(path, 'w') as f:
            f.write('objname')
            for category in catagories:
                f.write(delimiter+category)
            f.write('\n')
            for objname in data_dict:
                f.write(objname)
                for category in catagories:
                    f.write(delimiter+str(data_dict[objname][category]))
                f.write('\n')
        print('[+++] Files packed to', path)
        return
    else:
        print("[!!!] Invalid packing option! ['pack'/'unpack']")
        return
def data_proccesser(data_set, individual='', mag_unc_max=0, flux_unc_max=0, quiet=False):
    print('[+++] Processing '+data_set+' data...')

    # Getting constants and paths
    const = get_constants()

    # Check quiet
    if quiet:
        sys.stdout = open(os.devnull,
                          'w')

    objs = {}
    if data_set == 'CSP':
        files = glob.glob(const['csp_data_loc'] + '*.txt')
        if len(individual) > 0 and os.path.isfile(individual):
            files = [individual]

        csp_sn1a91bg_only = True
        names91bglike = np.genfromtxt(const['csp91bglike_txt'], delimiter=', ', dtype=str)
        for file in files:
            with open(file, 'r') as f:
                hdr = f.readline()[:-1].split(' ')
                objname, z, ra, dec = hdr[0][2:], hdr[1], hdr[2], hdr[3]

                if csp_sn1a91bg_only and (objname not in names91bglike):
                    continue

                print('---------------------------------------------------------------------------------------------------')
                print('[', files.index(file) + 1, '/', len(files), ']')
                print('[' + str(ra) + ', ' + str(dec) + ']', '--', objname, '|', z)

                objs.update({objname: {'ra': float(ra), 'dec': float(dec), 'z': float(z),
                                       'zp': np.array([]), 'filters': np.array([]), 'time': np.array([]),
                                       'flux': np.array([]), 'dflux': np.array([]), 'mag': np.array([]), 'dmag': np.array([])}})
                filter = ''
                for line in f.readlines():
                    data_line = line[:-1].split(' ')
                    if len(data_line) == 2:
                        filter = data_line[1]
                        continue
                    if len(data_line) == 3:
                        m, m_err = float(data_line[1]), float(data_line[2])
                        zp = float(const['csp_zpts_'+filter])
                        flux = 10**((m - zp) / 2.5) # dropping the negative makes it look more accurate, idk y
                        flux_err = 10**(m_err / 2.5)

                        objs[objname]['zp'] = np.append(objs[objname]['zp'], zp)
                        objs[objname]['filters'] = np.append(objs[objname]['filters'], filter)
                        objs[objname]['time'] = np.append(objs[objname]['time'], data_line[0])
                        objs[objname]['flux'] = np.append(objs[objname]['flux'], flux)
                        objs[objname]['dflux'] = np.append(objs[objname]['dflux'], flux_err)
                        objs[objname]['mag'] = np.append(objs[objname]['mag'], data_line[1])
                        objs[objname]['dmag'] = np.append(objs[objname]['dmag'], data_line[2])
    elif data_set == 'ATLAS':
        files = glob.glob(const['atlas_data_loc']+'*.txt')
        if len(individual) > 0 and os.path.isfile(individual):
            files = [individual]

        for file in files:
            print('-----------------------------------------------------------------------------------------------')
            print('[', files.index(file) + 1, '/', len(files), ']')

            # Pull data
            data = np.genfromtxt(file, delimiter=',', dtype=str, skip_header=1)
            if len(data) == 0:
                print('File empty!')
                continue

            # Get details
            ra, dec = np.average(data[:, 1].astype(float)), np.average(data[:, 2].astype(float))
            objname, z = TNS_objname_z(ra, dec)
            if z == 'None':
                z = np.nan
            else:
                z = float(z)
            print('['+str(ra)+', '+str(dec)+']', '--', objname, '|', z)

            # Update dictionary
            objs.update({objname: {'ra': ra, 'dec': dec, 'z': z, 'zp': data[:, 7].astype(float),
                                   'filters': data[:, 6], 'time': data[:, 8],
                                   'flux': data[:, 16], 'dflux': data[:, 17],
                                   'mag': data[:, 3], 'dmag': data[:, 4]}})
    elif data_set == 'ZTF':
        files = glob.glob(const['ztf_data_loc']+'*.txt')
        if len(individual) > 0 and os.path.isfile(individual):
            files = [individual]

        for file in files:
            print('-----------------------------------------------------------------------------------------------')
            print('[', files.index(file) + 1, '/', len(files), ']')

            # Pull data
            data = np.genfromtxt(file, delimiter=None, dtype=str, skip_header=56)
            if len(data) == 0:
                print('File empty!')
                continue

            # Header
            with open(file, 'r') as f:
                for i in range(55):
                    if i == 3:
                        ra = float(f.readline().split(' ')[-2])
                    elif i == 4:
                        dec = float(f.readline().split(' ')[-2])
                    else:
                        f.readline()
                hdr = f.readline()[1:-1].split(', ')

            # Get details
            objname, z = TNS_objname_z(ra, dec)
            print('['+str(ra)+', '+str(dec)+']', '--', objname, '|', z)

            # Get magnitudes m = -2.5log(F) + zp
            mags, dmags = np.array([]), np.array([])
            for i in range(len(data[:, 24])):
                n_flux, n_dflux, n_zp = data[:, 24][i], data[:, 25][i], data[:, 20].astype(float)[i]
                if n_flux == 'null' or n_dflux == 'null':
                    mags = np.append(mags, np.nan)
                    dmags = np.append(dmags, np.nan)
                else:
                    m = (-2.5*np.log10(float(n_flux))) + n_zp
                    m_err = 2.5*np.log10(float(n_dflux))
                    mags = np.append(mags, m)
                    dmags = np.append(dmags, m_err)

            # Update dictionary
            objs.update({objname: {'ra': ra, 'dec': dec, 'z': float(z), 'zp': data[:, 20].astype(float),
                                   'filters': data[:, 4], 'time': data[:, 22],
                                   'flux': data[:, 24], 'dflux': data[:, 25],
                                   'mag': mags, 'dmag': dmags}})
    else:
        raise ValueError("Data set not recognized ['CSP'/'ATLAS'/'ZTF']")
        return

    # Clean Data -- ['zp', 'filter', 'time', 'flux', 'dflux', 'mag', 'dmag']
    for obj in list(objs.keys()):
        new_zp, new_filter, new_time, new_mag, new_mag_unc, new_flux, new_flux_unc = (
            np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]))
        for n in range(len(objs[obj]['mag'])):
            n_zp, n_filter, n_time, n_mag, n_mag_unc, n_flux, n_flux_unc = (
                objs[obj]['zp'][n], objs[obj]['filters'][n], objs[obj]['time'][n],
                objs[obj]['mag'][n], objs[obj]['dmag'][n], objs[obj]['flux'][n], objs[obj]['dflux'][n])

            # Normallize Data
            norm_failed = False
            n_mag = str(n_mag).replace('>', '')
            for n in [n_mag, n_mag_unc, n_flux, n_flux_unc]:
                if str(n) == 'None' or str(n) == 'nan' or str(n) == 'null':
                    norm_failed = True
            if norm_failed:
                continue

            # Clean Data
            clean_failed = False
            for n in [n_mag, n_mag_unc, n_flux, n_flux_unc]:
                if float(n) <= 0:
                    clean_failed = True
            if clean_failed:
                continue
            if (mag_unc_max != 0) and (float(n_mag_unc) > mag_unc_max):
                continue
            if (flux_unc_max != 0) and (float(n_flux_unc) > flux_unc_max):
                continue

            # print('flux:', n_flux, '\t| dflux:', n_flux_unc, '\t| mag:', n_mag, '\t| dmag:', n_mag_unc)

            new_zp = np.append(new_zp, float(n_zp))
            new_filter = np.append(new_filter, n_filter)
            new_time = np.append(new_time, float(n_time))
            new_mag = np.append(new_mag, float(n_mag))
            new_mag_unc = np.append(new_mag_unc, float(n_mag_unc))
            new_flux = np.append(new_flux, float(n_flux))
            new_flux_unc = np.append(new_flux_unc, float(n_flux_unc))

        # Update object
        (objs[obj]['zp'], objs[obj]['filters'], objs[obj]['time'], objs[obj]['mag'], objs[obj]['dmag'],
         objs[obj]['flux'], objs[obj]['dflux']) = (
            new_zp, new_filter, new_time, new_mag, new_mag_unc, new_flux, new_flux_unc)

    # Restore print statements
    sys.stdout = sys.__stdout__

    return objs
def host_mass(dict_path, save_loc='../default/', keep_data=True, update_saved=False, use_mass_key=True):
    data = dict_handler(choice='unpack', path=dict_path)
    all_mass, all_mass_err = [], []
    GHOST_DATA = get_constants()['ghost_data_loc']
    print('[+++] Finding host galaxy mass using GHOST...')

    # Get mass key
    mass_key = {}
    if use_mass_key:
        with open(get_constants()['mass_key_txt'], 'r') as f:
            temp = f.readlines()
            for line in temp:
                line = line[:-1].split(', ')
                if len(line) != 3:
                    continue
                mass_key.update({line[0]: {'mass': line[1], 'mass_err': line[2]}})

    failed_host_masses, err_messages = [], []
    for obj in data:
        objname = obj
        if obj[:2] == 'SN':
            objname = obj[2:]

        print('\n[', list(data).index(obj) + 1, '/', len(data), ']', objname, '|', data[obj]['ra'], ',', data[obj]['dec'])
        print('-------------------------------------------------------------------------------------------------------')

        if objname in mass_key:
            print('Found object in mass key! Pulling...')
            all_mass.append(mass_key[objname]['mass'])
            all_mass_err.append(mass_key[objname]['mass_err'])
            if update_saved:
                data[obj].update({'host_mass': mass_key[objname]['mass'], 'host_mass_err': mass_key[objname]['mass_err']})
            continue

        ra, dec, z = float(data[obj]['ra']), float(data[obj]['dec']), float(data[obj]['z'])

        transient_position = SkyCoord(ra, dec, unit=u.deg)
        try:
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                host_data = getTransientHosts(transientCoord=[transient_position], transientName=[obj], verbose=False,
                                              starcut="gentle", savepath=save_loc + 'ghost_stuff/',
                                              GHOSTpath=GHOST_DATA)

            # If GLADE can't get magnitudes -- using NED name to get SDSS data
            print('Identified Host Galaxy:', host_data.loc[0, 'NED_name'])
            if np.isnan(host_data.loc[0, 'gKronMag']) or np.isnan(host_data.loc[0, 'iKronMag']):
                print('GLADE does not contain the g-mag & i-mag, searching SDSS...')
                host_position = SkyCoord(host_data['raMean'], host_data['decMean'], unit=u.deg)
                result = SDSS.query_crossid(host_position, photoobj_fields=['modelMag_g', 'modelMagErr_g', 'modelMag_i',
                                                                            'modelMagErr_i'])

                # If GLADE or SDSS has magnitudes -- using PanSTARRS
                if result == None:
                    print('GLADE & SDSS do not contain the g-mag & i-mag, searching PanSTARRS...')
                    catalog_data = Catalogs.query_region(transient_position, radius=2 * u.arcsec, catalog="Panstarrs")
                    # print(catalog_data['rMeanKronMag'], catalog_data['iMeanKronMag'])

                    if len(catalog_data) == 0 or isinstance(catalog_data['gMeanKronMagErr'].value[0],
                                                            np.ma.core.MaskedConstant) or isinstance(
                            catalog_data['iMeanKronMagErr'].value[0], np.ma.core.MaskedConstant):
                        print('[!!!!!] GLADE, SDSS, and PanSTARRS do not contain the g-mag & i-mag data...')
                        err_messages.append('GLADE, SDSS, and PanSTARRS do not contain the g-mag & i-mag data.')
                        gMag, iMag, iAbsMag = np.nan, np.nan, np.nan
                        gMagErr, iMagErr, iAbsMagErr = np.nan, np.nan, np.nan
                    else:
                        gMag, iMag, iAbsMag = catalog_data['gMeanKronMag'].value[0], catalog_data['gMeanKronMag'].value[
                            0], (catalog_data['iMeanKronMag'].value[0] - CURRENT_COSMO.distmod(z).value)
                        gMagErr, iMagErr, iAbsMagErr = catalog_data['gMeanKronMagErr'].value[0], \
                        catalog_data['iMeanKronMagErr'].value[0], catalog_data['iMeanKronMagErr'].value[0]
                else:
                    # SDSS Results
                    gMag, iMag, iAbsMag = result['modelMag_g'].value[0], result['modelMag_i'].value[0], (
                                result['modelMag_i'].value[0] - CURRENT_COSMO.distmod(z).value)
                    gMagErr, iMagErr, iAbsMagErr = result['modelMagErr_g'].value[0], result['modelMagErr_i'].value[0], \
                    result['modelMagErr_i'].value[0]
            else:
                # GLADE results
                gMag, iMag, iAbsMag = host_data['gKronMag'].loc[0], host_data['iKronMag'].loc[0], (
                            host_data['iKronMag'].loc[0] - CURRENT_COSMO.distmod(z).value)
                gMagErr, iMagErr, iAbsMagErr = host_data['gKronMagErr'].loc[0], host_data['iKronMagErr'].loc[0], \
                host_data['iKronMagErr'].loc[0]

            #  Mass Calculation -- Taylor et. al. 2011 -- eq. 8
            host_mass = (1.15 + (0.7 * (gMag - iMag)) - (0.4 * (iAbsMag)))

            # Error Propogation
            giMagErr = np.sqrt((gMagErr ** 2) + (iMagErr ** 2))
            host_mass_err = np.sqrt(((0.7 ** 2) * (giMagErr ** 2)) + ((0.4 ** 2) * (iAbsMagErr ** 2)))

        except Exception as error:
            err_messages.append('Target does not exsist in GLADE, SDSS, or PanSTARRS.')
            print(obj + ':', error)
            host_mass, host_mass_err = 0.00, 0.00

        if np.isnan(host_mass) == False and host_mass > 0.00:
            print('Success!', obj, 'host galaxy has a mass of:', host_mass, '+/-', host_mass_err, 'logM_* / [Msun]')
            all_mass.append(host_mass)
            all_mass_err.append(all_mass_err)
            # Update mass key
            if use_mass_key:
                with open(get_constants()['mass_key_txt'], 'a') as f:
                    print('Updating mass key with '+obj+'...')
                    f.write(objname+', '+str(host_mass)+', '+str(host_mass_err)+'\n')
            if update_saved:
                data[obj].update({'host_mass': host_mass, 'host_mass_err': host_mass_err})
        else:
            print('[!!!] Failed to find host galaxy!')
            if update_saved:
                data[obj].update({'host_mass': 0.00, 'host_mass_err': 0.00})
            failed_host_masses.append([obj, ra, dec, z])

    print('\nSuccessfully found mass of', len(all_mass), '/', len(data), 'host galaxies!')
    if len(failed_host_masses) > 0:
        print('************************************************************************************************')
        print('[!!!] Failed host mass calculations:')
        print('### [objname, ra, dec, z] ###')
        for fail in failed_host_masses:
            print(fail, '--', err_messages[failed_host_masses.index(fail)])
        print('************************************************************************************************')


    if not keep_data:
        print('Removing GHOST data...')
        shutil.rmtree(save_loc + 'ghost_stuff/')  # Clear messy data
        os.mkdir(save_loc + 'ghost_stuff/')
    if update_saved:
        print('Saving data to' + dict_path + '...')
        dict_handler(choice='pack', data_dict=data, path=dict_path)
    return all_mass, all_mass_err
def get_reviewed_fits(path=get_constants()['reviewed_fits_txt']):
    reviewed_fits = {'good': [], 'okay': [], 'bad': []}
    with open(path, 'r') as f:
        # reviewed_good_fits, reviewed_okay_fits, reviewed_bad_fits = [], [], []
        for container in [reviewed_fits['good'], reviewed_fits['okay'], reviewed_fits['bad']]:
            for i in f.readline().split(', ')[1:]:
                container.append(i[1:-1])
            container[-1] = container[-1][:-1]
    return reviewed_fits
def mass_step_calc(path, cut=10, ignore_fits=[], quiet=False):
    objs = dict_handler(choice='unpack', path=path)
    lower_mass, upper_mass, all_resid = np.array([]), np.array([]), np.array([])
    lower_mass_err, upper_mass_err, all_resid_err = np.array([]), np.array([]), np.array([])
    reviewed_fits = get_reviewed_fits()

    # Check quiet
    if quiet:
        sys.stdout = open(os.devnull, 'w')

    for obj in objs:
        if ('bad' in ignore_fits) and (obj in reviewed_fits['bad']):
            continue
        elif ('okay' in ignore_fits) and (obj in reviewed_fits['okay']):
            continue
        elif ('good' in ignore_fits) and (obj in reviewed_fits['good']):
            continue

        m = float(objs[obj]['host_mass'])
        m_err = float(objs[obj]['host_mass_err'])
        resid = float(objs[obj]['mu']) - CURRENT_COSMO.distmod(float(objs[obj]['z'])).value
        resid_err = float(objs[obj]['mu_err'])

        all_resid = np.append(all_resid, resid)
        all_resid_err = np.append(all_resid_err, resid_err)
        if m < cut:
            lower_mass = np.append(lower_mass, m)
            lower_mass_err = np.append(lower_mass_err, m_err)
        else:
            upper_mass = np.append(upper_mass, m)
            upper_mass_err = np.append(upper_mass_err, m_err)
    mass_step, mass_step_err = np.average(upper_mass) - np.average(lower_mass), abs(np.average(upper_mass_err) - np.average(lower_mass_err))

    print('***********************************************************************************************************')
    print('Mass Step (M_>'+str(cut)+' - M_<'+str(cut)+'):', mass_step, '+/-', mass_step_err)
    print('\t----------------------')
    print('\tLog Mass > '+str(cut)+':', round(np.average(upper_mass), 4), '+/-', abs(round(np.average(upper_mass_err), 4)))
    print('\tLog Mass < '+str(cut)+':', round(np.average(lower_mass), 4), '+/-', abs(round(np.average(lower_mass_err), 4)))
    print('\t----------------------')
    print('\tScatter:', round(np.std(all_resid), 4))
    print('\t----------------------')
    print('\tNumber of Targets:', len(all_resid))
    print('***********************************************************************************************************')

    # Restore print statements
    sys.stdout = sys.__stdout__

    return mass_step, mass_step_err
# ==================================================================================================================== #
def lc_plot(objs, y_type = 'flux', pause_time=2, color_wheel = ['orange', 'cyan', 'violet', 'red', 'blue'],
            quiet=False, save_plots=True, save_loc='../snpy/misc_plots/'):
    print('[+++] Plotting LC data...')
    color_wheel = [None, None, None, None, None, None, None, None, None, None, None]

    # Check quiet
    if quiet:
        sys.stdout = open(os.devnull, 'w')

    for obj in objs:
        print('-------------------------------------------------------------------------------------------------------')
        print('[', list(objs).index(obj)+1, '/', len(objs), '] -', obj)

        plt.figure(figsize=(12, 6))

        # Get list of all filters
        filter_wheel = []
        for f_w in objs[obj]['filters']:
            if f_w not in filter_wheel:
                filter_wheel.append(f_w)

        for f_w in filter_wheel:
            f_indexs = np.where(objs[obj]['filters'] == f_w)[0]
            plt.errorbar(objs[obj]['time'][f_indexs], objs[obj][y_type][f_indexs], yerr=objs[obj]['d'+y_type][f_indexs],
                        color=color_wheel[filter_wheel.index(f_w)], label=f_w, fmt='o')
        systime.sleep(pause_time)

        if y_type == 'mag':
            plt.gca().invert_yaxis()
        plt.title(obj)
        plt.xlabel('JD'); plt.ylabel(y_type)
        plt.legend()
        if save_plots:
            plt.savefig(save_loc + obj + '_lc.png')
            print(obj, '-- Plot saved to', save_loc + obj + '_lc.png')
        plt.show()
        plt.close()


    # Restore print statements
    sys.stdout = sys.__stdout__

    return
def lc_replot(lc_path, save_plot=False, save_loc='../default/', colors=None, spread=None, stacked=True):
    n_s = snpy.get_sn(lc_path)
    # if colors is None:
    #     plt_args = {'single': stacked, 'colors': colors}
    # else:
    plt_args = {'single': stacked}

    if spread != None:
        plt_args.update({'xrange': (n_s.parameters['Tmax']-spread[0], n_s.parameters['Tmax']+spread[1])})
    if save_plot:
        plt_args.update({'outfile': save_loc + n_s.name + '_snpylc.png'})
    n_s.plot(**plt_args)
    plt.show()
    return
def residual_plotter(path, x_params, sigma=[1,1], labels=False, raw=False, extra_info=False, ignore_type=[], save_plot=False, save_loc='../default/'):
    # Get reviewed fits
    with open(get_constants()['reviewed_fits_txt'], 'r') as f:
        reviewed_good_fits, reviewed_okay_fits, reviewed_bad_fits = [], [], []
        for container in [reviewed_good_fits, reviewed_okay_fits, reviewed_bad_fits]:
            for i in f.readline().split(', ')[1:]:
                container.append(i[1:-1])
            container[-1] = container[-1][:-1]

    # Pull data from saved text
    objs = dict_handler(choice='unpack', path=path)

    fig, axs = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'width_ratios': [10, 1]}, constrained_layout=True)
    atlas_key, ztf_key, csp_key = True, True, True
    colors = {'ATLAS': 'blue', 'ZTF': 'red', 'CSP':'orange'}
    mu_res_hist = []
    for obj in objs:
        if x_params[0] == 'host_mass':
            n_x_err = float(objs[obj][x_params[0]+'_err'])
        else:
            n_x_err = 0.00
        n_x = float(objs[obj][x_params[0]])
        n_mu_res = float(objs[obj]['mu']) - CURRENT_COSMO.distmod(float(objs[obj]['z'])).value
        n_mu_err = float(objs[obj]['mu_err'])

        # Clean data
        if not raw:
            objname = obj
            if objname[:2] == 'SN':
                objname = objname[2:]
            if ('good' in ignore_type) and (objname in reviewed_good_fits):
                continue
            if ('okay' in ignore_type) and (objname in reviewed_okay_fits):
                continue
            if ('bad' in ignore_type) and (objname in reviewed_bad_fits):
                continue
            if x_params[0] == 'host_mass' and n_x == 0:  # Remove null masses
                continue

        # Save for histogram
        mu_res_hist.append(n_mu_res)

        # Plot points
        if path.split('/')[-1][:-10] == 'combined':
            if atlas_key and (objs[obj]['origin'] == 'ATLAS'):
                axs[0].errorbar(n_x, n_mu_res, yerr=n_mu_err * sigma[0], xerr=n_x_err * sigma[1],
                                fmt='o', color=colors[objs[obj]['origin']], label=objs[obj]['origin'])
                atlas_key = False
            elif ztf_key and (objs[obj]['origin'] == 'ZTF'):
                axs[0].errorbar(n_x, n_mu_res, yerr=n_mu_err * sigma[0], xerr=n_x_err * sigma[1],
                                fmt='o', color=colors[objs[obj]['origin']], label=objs[obj]['origin'])
                ztf_key = False
            elif csp_key and (objs[obj]['origin'] == 'CSP'):
                axs[0].errorbar(n_x, n_mu_res, yerr=n_mu_err * sigma[0], xerr=n_x_err * sigma[1],
                                fmt='o', color=colors[objs[obj]['origin']], label=objs[obj]['origin'])
                csp_key = False
            else:
                axs[0].errorbar(n_x, n_mu_res, yerr=n_mu_err * sigma[0], xerr=n_x_err * sigma[1],
                                fmt='o', color=colors[objs[obj]['origin']])
        else:
            axs[0].errorbar(n_x, n_mu_res, yerr=n_mu_err * sigma[0], xerr=n_x_err * sigma[1],
                            fmt='o')

        # Labels
        if labels:
            axs[0].text(n_x, n_mu_res, obj, size='x-small', va='top')

    # Histogram
    axs[1].hist(mu_res_hist, bins=40, orientation="horizontal")

    # Formatting
    ylimiter = np.max(np.abs(mu_res_hist))+0.5
    axs[0].set_ylim(-ylimiter, ylimiter); axs[1].set_ylim(-ylimiter, ylimiter)

    if extra_info:
        fig.suptitle("Hubble Residuals vs. " + x_params[1] + " of '"+(path.split('/')[-1].split('_')[0]).upper()+"' 91bg-like SNe Ia\n" +  # Figure Title
                     'Dist. Sigma: ' + str(sigma[0]) + ' | ' + x_params[1] + ' Sigma: ' + str(sigma[0]) +
                     ' | Scatter: ' + str(round(np.std(mu_res_hist), 2)) + ' | # of pts: ' + str(len(mu_res_hist)), size='medium')
    else:
        fig.suptitle("Hubble Residuals vs. " + x_params[1] + " of '"+(path.split('/')[-1][:-10]).upper()+"' 91bg-like SNe Ia")
    axs[0].set(xlabel=x_params[1], ylabel='Hubble Residuals')  # Sub-plot Labels
    axs[1].get_yaxis().set_visible(False) # Turn off y-axis labels
    if path.split('/')[-1][:-10] == 'combined':
        axs[0].legend()
    if save_plot:
        plt.savefig(save_loc+'hubble_res_v_'+x_params[1]+'.png')
    plt.show()
    return
def snpy_histogram_plotter(path, raw=False, save_plot=False, save_loc='../default/', ignore_type=[], param_bins=[None, None, None, None, None]):
    # Get reviewed fits
    with open(get_constants()['reviewed_fits_txt'], 'r') as f:
        reviewed_good_fits, reviewed_okay_fits, reviewed_bad_fits = [], [], []
        for container in [reviewed_good_fits, reviewed_okay_fits, reviewed_bad_fits]:
            for i in f.readline().split(', ')[1:]:
                container.append(i[1:-1])
            container[-1] = container[-1][:-1]

    # Pull data
    objs = dict_handler(path=path, choice='unpack')
    mu, st, Tmax, EBVhost, host_mass = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
    for obj in objs:
        # Clean data
        if not raw:
            objname = obj
            if objname[:2] == 'SN':
                objname = objname[2:]
            if ('good' in ignore_type) and (objname in reviewed_good_fits):
                continue
            if ('okay' in ignore_type) and (objname in reviewed_okay_fits):
                continue
            if ('bad' in ignore_type) and (objname in reviewed_bad_fits):
                continue

        mu = np.append(mu, float(objs[obj]['mu']))
        st = np.append(st, float(objs[obj]['st']))
        Tmax = np.append(Tmax, float(objs[obj]['Tmax']))
        EBVhost = np.append(EBVhost, float(objs[obj]['EBVhost']))
        host_mass = np.append(host_mass, float(objs[obj]['host_mass']))

    # Plot
    fig, ax = plt.subplots(1, 5, figsize=(16, 4), layout='constrained')
    params = [mu, st, Tmax, EBVhost, host_mass]
    param_names = ['mu', 'st', 'Tmax', 'EBVhost', 'host_mass']
    param_bins = [45, 45, 45, 45, 45]
    for i in range(len(params)):
        ax[i].hist(params[i], bins=param_bins[i])
        if i != 0:
            ax[i].get_yaxis().set_visible(False)
        ax[i].set_xlabel(param_names[i])

    plt.suptitle("Parameters for '" + path.split('/')[-1].split('_')[0].upper()
                 + "' data\n Number of Transients: " + str(len(objs)), fontsize=20)
    if save_plot:
        print('Saved figure to... ', save_loc+path.split('/')[-1].split('_')[0]+'_hist.png')
        plt.savefig(save_loc+path.split('/')[-1].split('_')[0]+'_hist.png')
    plt.show()

    return