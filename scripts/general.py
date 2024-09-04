import os
import sys
import glob
import shutil
import time as systime
import numpy as np

# from astro_ghost.ghostHelperFunctions import getTransientHosts
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u

from astroquery.sdss import SDSS
from astroquery.mast import Catalogs
import warnings

def current_cosmo(H0=70, O_m=0.3):
    return FlatLambdaCDM(H0, O_m)
def get_constants(cosntant_loc='txts/constants.txt'):
    CONSTANTS = {}
    with open(cosntant_loc, 'r') as f:
        temp = f.readlines()
        for line in temp:
            line = line[:-1].split(', ')
            if len(line) != 2:
                continue
            CONSTANTS.update({line[0]: line[1]})
    return CONSTANTS
def get_APIkeys(apikeys_loc='txts/APIkeys.txt'):
    if not os.path.isfile(apikeys_loc):
        raise FileNotFoundError('[!!!] API keys file not found!')
    APIKEY = {}
    with open(apikeys_loc, 'r') as f:
        temp = f.readlines()
        for line in temp:
            line = line[:-1].split(', ')
            if len(line) != 2:
                continue
            APIKEY.update({line[0]: line[1]})
    return APIKEY
def TNS_objname_z(obj_ra, obj_dec, attempts=10, use_keys=True):
    from scripts import tnsAPI
    print('-----------------------------------------------------------------------------------------------------------')
    print('Searching TNS for ['+str(obj_ra)+', '+str(obj_dec)+']...')
    with open('txts/APIkeys.txt', 'r') as f:
        tns_bot_id = f.readline().split(', ')[-1][:-1]
        tns_bot_name = f.readline().split(', ')[-1][:-1]
        tns_bot_api_key = f.readline().split(', ')[-1][:-1]
    tns_key = dict_handler(choice='unpack', path=get_constants()['tns_key_txt'])
    obj_name, obj_z = '', 0.00

    # Check key
    print('Checking TNS key...')
    if use_keys and (str(obj_ra) in tns_key):
        obj_name, obj_z = tns_key[str(obj_ra)]['objname'], tns_key[str(obj_ra)]['z']
    else:
        # Query TNS key
        print('Not in TNS key, querying TNS...')
        try:
            # Code abridged from David's code
            tns_bot_id, tns_bot_name, tns_bot_api_key = '73181', 'YSE_Bot1', '0d771345fa6b876a5bb99cd5042ab8b5ae91fc67'
            headers = tnsAPI.build_tns_header(tns_bot_id, tns_bot_name)
            tns_api_url = f"https://www.wis-tns.org/api/get"
            search_tns_url = tnsAPI.build_tns_url(tns_api_url, mode="search")
            get_tns_url = tnsAPI.build_tns_url(tns_api_url, mode="get")
            search_data = tnsAPI.build_tns_search_query_data(tns_bot_api_key, obj_ra, obj_dec)
            transients = tnsAPI.rate_limit_query_tns(search_data, headers, search_tns_url)
            get_data = tnsAPI.build_tns_get_query_data(tns_bot_api_key, transients[0])
            details = tnsAPI.rate_limit_query_tns(get_data, headers, get_tns_url)
            obj_name, obj_z = details['objname'], details['redshift']
            with open(get_constants()['tns_key_txt'], 'a') as f:
                f.write(str(obj_ra) + ', ' + str(obj_dec) + ', ' + obj_name + ', ' + str(obj_z) + '\n') # Save to key
        except Exception as error:
            if attempts <= 0:
                raise RuntimeError('TNS completly timed out.')
            print(f'{attempts} attempts remaining')
            systime.sleep(5)
            obj_name, obj_z = TNS_objname_z(obj_ra, obj_dec, attempts=attempts-1, use_keys=use_keys)

    return obj_name, obj_z
def dict_handler(data_dict={}, choice='', path='default/dict.txt', delimiter=', '):
    if choice == 'unpack':
        print('[+++] Unpacking objects from '+path+'...')
        with open(path, 'r') as f:
            hdr = f.readline()[:-1].split(delimiter)
        data = np.genfromtxt(path, delimiter=delimiter, dtype=str, skip_header=1)
        if len(data) == 0:
            return {}
        temp_objs = {}

        if len(np.shape(data)) == 1:
            iter_period = len(data[0])
        else:
            iter_period = len(data[:, 0])

        for i in range(iter_period):
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
    elif choice == 'arrays':
        print('[+++] Unpacking objects from '+path+' into arraies...')
        with open(path, 'r') as f:
            hdr = f.readline()[:-1].split(delimiter)
        data = np.genfromtxt(path, delimiter=delimiter, dtype=str, skip_header=1)
        if len(data) == 0:
            return {}
        temp_arrs = {}
        for i in range(len(hdr)):
            temp_arrs.update({hdr[i]: data[:, i]})
        return temp_arrs
    else:
        print("[!!!] Invalid packing option! ['arrays'/'pack'/'unpack']")
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
            with warnings.catch_warnings():
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
        print(individual)
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
def host_mass(dict_path, save_loc='default/', keep_data=True, update_saved=False, use_mass_key=True):
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
                            0], (catalog_data['iMeanKronMag'].value[0] - current_cosmo().distmod(z).value)
                        gMagErr, iMagErr, iAbsMagErr = catalog_data['gMeanKronMagErr'].value[0], \
                        catalog_data['iMeanKronMagErr'].value[0], catalog_data['iMeanKronMagErr'].value[0]
                else:
                    # SDSS Results
                    gMag, iMag, iAbsMag = result['modelMag_g'].value[0], result['modelMag_i'].value[0], (
                                result['modelMag_i'].value[0] - current_cosmo().distmod(z).value)
                    gMagErr, iMagErr, iAbsMagErr = result['modelMagErr_g'].value[0], result['modelMagErr_i'].value[0], \
                    result['modelMagErr_i'].value[0]
            else:
                # GLADE results
                gMag, iMag, iAbsMag = host_data['gKronMag'].loc[0], host_data['iKronMag'].loc[0], (
                            host_data['iKronMag'].loc[0] - current_cosmo().distmod(z).value)
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
        resid = float(objs[obj]['mu']) - current_cosmo().distmod(float(objs[obj]['z'])).value
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
