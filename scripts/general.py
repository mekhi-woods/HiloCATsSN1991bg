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
def TNS_unpack(obj_ra, obj_dec, check_key=True):
    print('Searching TNS for ['+str(obj_ra)+', '+str(obj_dec)+']...')
    tns_key = np.genfromtxt(get_constants()['tns_key_txt'], dtype=str, delimiter=', ', skip_header=1)

    # Check TNS key
    name_found = False
    if len(tns_key) > 0:
        print('Checking TNS key...')
        objnames, ra, dec, z = tns_key[:, 2], tns_key[:, 0].astype(float), tns_key[:, 1].astype(float), tns_key[:, 3]
        for n in range(len(objnames)):
            if abs(obj_ra - ra[n]) < 0.01:
                obj_name, obj_z = objnames[n], z[n]
                print('Found! ['+str(obj_ra)+', '+str(obj_dec)+'] -- ', obj_name, '|', obj_z)
                name_found = True
                break

    # If no check or not found, query TNS
    if name_found == False:
        print('Querying TNS...')
        details = TNS_details(obj_ra, obj_dec)
        obj_name, obj_z = details['objname'], details['redshift']
        print('Found! ['+str(obj_ra)+', '+str(obj_dec)+'] -- ', obj_name, ',', obj_z)
        # Save to key
        with open(get_constants()['tns_key_txt'], 'a') as f:
            f.write(str(obj_ra) + ', ' + str(obj_dec) + ', ' + obj_name + ', ' + str(obj_z) + '\n')

    return obj_name, obj_z
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
        print('[+++] Unpacking objects from text file...')
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
def handle_problem_children(state, problem_c=None):
    path = get_constants()['prob_child_txt']

    if state == 'READ':
        # Read problem children
        problem_c = np.genfromtxt(path, dtype=str)
        return problem_c
    elif state == 'WRITE':
        # Write problem children
        problem_c = np.unique(problem_c)
        with open(path, 'w') as f:
            for c in problem_c:
                f.write(c+'\n')
        return None
    else:
        raise Exception("Invalid state: '"+state+"' [READ/WRITE]")
# ==================================================================================================================== #
def data_proccesser(data_set = 'ZTF', mag_unc_max=1, quiet=False):
    print('[+++] Processing '+data_set+' data...')

    # Getting constants and paths
    const = get_constants()

    # Check quiet
    if quiet:
        sys.stdout = open(os.devnull, 'w')

    # Chose data set
    if data_set == 'ATLAS':
        hdr_len = 1
        data_delim, hdr_delim = ',', ','
        files = glob.glob(const['atlas_data_loc']+'*.txt')
        d_i = {'mag': 3, 'dmag': 4, 'filter': 6, 'zp': 7, 'time': 8, 'flux': 16, 'dflux': 17, 'uJy': 24, 'duJy': 25}
    elif data_set == 'ZTF':
        hdr_len = 56
        data_delim, hdr_delim = None, ', '
        files = glob.glob(const['ztf_data_loc']+'*.txt')
        d_i = {'filter': 4, 'zp': 20, 'time': 22, 'flux': 24, 'dflux': 25}
    else:
        print('Data set not recognized')
        return None

    # Getting data & header
    data_objs = {}
    for file in files:
        print('-------------------------------------------------------------------------------------------------------')
        print('[', files.index(file)+1, '/', len(files), ']')

        # Retrieve Data
        data = np.genfromtxt(file, dtype=str, skip_header=hdr_len, delimiter=data_delim)

        # Check if file is empty
        if len(data) <= 0:
            print('[!!!] File empty!')
            continue

        # Header
        with open(file, 'r') as f:
            if hdr_len > 1:
                for i in range(hdr_len-1):
                    if data_set == 'ZTF' and i == 3:
                        obj_ra = float(f.readline().split(' ')[-2])
                    elif data_set == 'ZTF' and i == 4:
                        obj_dec = float(f.readline().split(' ')[-2])
                    else:
                        f.readline()
            hdr = f.readline()[1:-1].split(hdr_delim)

        # RA & Dec for ATLAS
        if data_set == 'ATLAS':
            obj_ra, obj_dec = np.average(data[:, 1].astype(float)), np.average(data[:, 2].astype(float))

        # Get name and redshift from TNS
        run = True
        obj_name = ''
        while run:
            try:
                obj_name, obj_z = TNS_unpack(obj_ra, obj_dec)
                if len(obj_name) > 0:
                    run = False
            except Exception as error:
                if error == NameError:
                    run = False
                print('***********************************************************************************************')
                print(error)
                print('TNS timed out, pausing for 5 seconds...')
                print('***********************************************************************************************')
                systime.sleep(5)
        print(obj_name, '|', file.split('/')[-1][:-4])
        print('File @', file)

        # Object Creation
        obj_filters, obj_time, obj_flux, obj_flux_unc, obj_mag, obj_mag_unc, obj_zp = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
        for i in range(len(data[:, 0])):
            # Flux to magnitude
            if data_set == 'ZTF':
                if data[i, d_i['flux']] == 'null':
                    continue
                else:
                    n_flux, n_flux_unc, n_flux_zpt = float(data[i, d_i['flux']]), float(data[i, d_i['dflux']]), float(data[i, 10]),
                    n_mag = -2.5 * np.log10(n_flux) + n_flux_zpt
                    n_mag_unc = (2.5 / (n_flux * np.log(10))) * n_flux_unc
                    if np.isnan(n_mag):
                        continue

            # Set variables
            if data_set == 'ATLAS':
                if data[i, d_i['mag']].split('>')[-1] == 'None':
                    continue
                else:
                    n_mag = float(data[i, d_i['mag']].split('>')[-1])
                    n_mag_unc = float(data[i, d_i['dmag']])
                n_flux, n_flux_unc = float(data[i, d_i['flux']]), float(data[i, d_i['dflux']])
                n_filter = data[i, d_i['filter']]
                n_time = float(data[i, d_i['time']])

            # Clean data
            if n_mag <= 0:  # Negatives mags
                continue
            if n_mag_unc <= 0:  # Negatives mag errors
                continue
            if n_flux <= 0: # Negatives fluxes
                continue
            if n_flux_unc <= 0:  # Negatives flux errors
                continue
            if (mag_unc_max != 0) and (n_mag_unc > mag_unc_max):
                continue

            # Commit variables
            obj_filters = np.append(obj_filters, data[i, d_i['filter']])
            obj_time = np.append(obj_time, float(data[i, d_i['time']]))
            obj_zp = np.append(obj_zp, float(data[i, d_i['zp']]))
            obj_flux = np.append(obj_flux, n_flux)
            obj_flux_unc = np.append(obj_flux_unc, n_flux_unc)
            obj_mag = np.append(obj_mag, n_mag)
            obj_mag_unc = np.append(obj_mag_unc, n_mag_unc)

        # Create Object
        obj = {'ra': obj_ra, 'dec': obj_dec, 'z': obj_z, 'filters': obj_filters, 'time': obj_time, 'zp': obj_zp,
               'flux': obj_flux, 'dflux': obj_flux_unc, 'mag': obj_mag, 'dmag': obj_mag_unc}
        print('Retrived: ' + obj_name + ' | ra: '+str(obj_ra)+' \ dec: '+str(obj_dec))

        # Commit Object
        data_objs.update({obj_name: obj})

    # Restore print statements
    sys.stdout = sys.__stdout__

    return data_objs
def alt_data_proccesser(data_set = 'ZTF', individual='', quiet=False):
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
                                   'flux': data[:, 16], 'dflux': data[:, 16],
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
    new_zp, new_filter, new_time, new_mag, new_mag_unc, new_flux, new_flux_unc = (
        np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]))
    for obj in list(objs.keys()):
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
def sample_cutter(save_loc):
    print('[+++] Cutting sample for ' + save_loc.split('_')[0].upper() + ' data...')

    path = glob.glob(get_constants()[save_loc] + '*_saved.txt')[0]
    objs = dict_handler(choice='unpack', path=path)
    i = 0
    new_objs = {}
    for obj in objs:
        print('[' + str(list(objs.keys()).index(obj) + 1) + '/' + str(len(objs)) + '] -- ' + obj)
        print('---------------------------------------------------------------------------------------------------')
        resid = float(objs[obj]['mu']) - CURRENT_COSMO.distmod(float(objs[obj]['z'])).value
        resid -= np.median(resid)

        if float(objs[obj]['EBVhost']) < -0.2 or float(objs[obj]['EBVhost']) > 0.2:
            print('[!!!] EBVhost out of range.')
            continue
        if float(objs[obj]['EBVhost_err']) > 0.1:
            print('[!!!] EBVhost errors out of range.')
            continue

        if float(objs[obj]['st']) < 0.3 or float(objs[obj]['st']) > 1.0:
            print('[!!!] Stretch out of range.')
            continue
        if float(objs[obj]['st_err']) > 0.1:
            print('[!!!] Stretch error out of range.')
            continue

        if float(objs[obj]['Tmax_err']) > 1:
            print('[!!!] Maximum time error out of range.')
            continue

        if float(objs[obj]['z']) < 0.015:
            print('[!!!] Redshift out of range.')
            continue
        i = i + 1

        # Save obj to new dict
        new_objs.update({obj: objs[obj]})

    dict_handler(data_dict=new_objs, choice='pack', path=get_constants()[save_loc] + save_loc.split('_')[0] + '_saved_cut.txt')
    return
# ==================================================================================================================== #
def write_ASCII(objs, filter_set, save_loc):
    print('[+++] Saving data to ASCII files for SNooPy...')
    for obj in objs:
        with open(save_loc + obj + '_snpy.txt', 'w') as f:
            # Line 1 -- Objname, Helio-Z, RA, Dec (Ex. SN1981D 0.005871 50.65992 -37.23272)
            f.write(str(obj)+' '+str(objs[obj]['z'])+' '+str(objs[obj]['ra'])+' '+str(objs[obj]['dec'])+'\n')
            for f_w in filter_set:
                f_indexs = np.where(objs[obj]['filters'] == f_w)[0]
                f.write('filter '+filter_set[f_w]+'\n')
                for i in f_indexs:
                    # filter photometry block -- Date (JD/MJD), mag, err (i.e. 674.8593 12.94 0.11)
                    f.write(str(objs[obj]['time'][i])+'\t'+str(objs[obj]['mag'][i])+'\t'+str(objs[obj]['dmag'][i])+'\n')
    return
def snpy_fit(paths, save_loc, use_saved=False, snpy_plots=True, save_plots=True, quiet=False):
    print('[+++] Fitting data with SNooPy...')

    # Check quiet
    if quiet:
        sys.stdout = open(os.devnull, 'w')

    objParams = {}
    for path in paths:
        objname = path.split('/')[-1][:-9]
        print('-------------------------------------------------------------------------------------------------------')
        print('[', paths.index(path)+1, '/', len(paths), '] --', objname)

        # Check saved fits
        saved_found = False
        if use_saved:
            for saved_fit in glob.glob(save_loc+'/models/*.snpy'):
                if objname == saved_fit.split('/')[-1][:-16]:
                    print('Saved fit found for '+objname+'...')
                    saved_found = True
                    n_s = snpy.get_sn(saved_fit)
                    mjds, mjde = 999999999, 0
                    for filter in list(n_s.data.keys()):
                        n_min, n_max = min(n_s.data[filter].MJD), max(n_s.data[filter].MJD)
                        if mjds > n_min:
                            mjds = n_min
                        if mjde < n_max:
                            mjde = n_max
                    objParams.update({objname: {'ra': n_s.ra, 'dec': n_s.decl, 'z': n_s.z, 'MJDs': mjds, 'MJDe': mjde,
                                                'mu': n_s.parameters['DM'], 'st': n_s.parameters['st'],
                                                'Tmax': n_s.parameters['Tmax'], 'EBVhost': n_s.parameters['EBVhost'],
                                                'mu_err': n_s.errors['DM'], 'st_err': n_s.errors['st'],
                                                'Tmax_err': n_s.errors['Tmax'], 'EBVhost_err': n_s.errors['EBVhost']}})
                    if snpy_plots:
                        n_s.plot()
                        plt.show()
                        systime.sleep(3)
                    print(objParams[objname])
        if saved_found:
            continue

        problem_children = handle_problem_children(state='READ') # Open problem children

        # Load data
        try:
            n_s = snpy.get_sn(path)
            n_s.choose_model('EBV_model2', stype='st')
            n_s.set_restbands()  # Auto pick appropriate rest-bands
        except:
            print('[!!!] Failed to load ASCII file')
            continue

        # Remove empty filters -- fix for 'ValueError: attempt to get argmin of an empty sequence'
        for filter in list(n_s.data.keys()):
            if len(n_s.data[filter].magnitude) == 0:
                del n_s.data[filter]
        print('Best filters:', list(n_s.data.keys()))

        run = True
        while run:
            try:
                # Minimum MJD and Max MJD
                mjds, mjde = 999999999, 0
                for filter in list(n_s.data.keys()):
                    n_min, n_max = min(n_s.data[filter].MJD), max(n_s.data[filter].MJD)
                    if mjds > n_min:
                        mjds = n_min
                    if mjde < n_max:
                        mjde = n_max

                n_s.fit(bands=None, dokcorr=True, k_stretch=False, reset_kcorrs=True, **{'mangle': 1, 'calibration': 0})
                n_s.save(save_loc + 'models/' + objname + '_EBV_model2.snpy')
                run = False

                if snpy_plots:
                    n_s.plot(outfile=save_loc + 'plots/' + objname + '_snpyplots.png')
                    plt.show()
                plt.close()

                objParams.update({objname : {'ra': n_s.ra, 'dec': n_s.decl, 'z': n_s.z, 'MJDs': mjds, 'MJDe': mjde,
                                             'mu': n_s.parameters['DM'], 'st': n_s.parameters['st'],
                                             'Tmax': n_s.parameters['Tmax'], 'EBVhost': n_s.parameters['EBVhost'],
                                             'mu_err': n_s.errors['DM'], 'st_err': n_s.errors['st'],
                                             'Tmax_err': n_s.errors['Tmax'], 'EBVhost_err': n_s.errors['EBVhost']}})
            except Exception as error:
                if 'All weights for filter' and 'are zero.' in str(error):
                    print('[!!!] Weights for filter', str(error).split(' ')[4], 'are zero. Removing...')
                    del n_s.data[str(error).split(' ')[4]]
                elif 'Error:  to solve for EBVhost, you need to fit more than one filter' in str(error):
                    print('[!!!]', error)
                    problem_children = np.append(problem_children, objname)
                    handle_problem_children(state='WRITE', problem_c=problem_children) # Commit problem children
                    run = False
                else:
                    print(error)
                    problem_children = np.append(problem_children, objname)
                    handle_problem_children(state='WRITE', problem_c=problem_children) # Commit problem children
                    run = False

    print('Successfully fit [', len(objParams), '/', len(paths), '] !')

    # Restore print statements
    sys.stdout = sys.__stdout__

    return objParams
def old_host_mass(dict_path, save_loc='../default/', keep_data=True, update_saved=False):
    cosmo = FlatLambdaCDM(70, 0.3) # Hubble Constant, Omega-Matter
    data = dict_handler(choice='unpack', path=dict_path)
    all_z, all_logstellarmass = [], []
    GHOST_DATA = get_constants()['ghost_data_loc']
    print('[+++] Finding host galaxy mass using GHOST...')

    i = 0
    for obj in data:

        # i += 1
        # if i < 1:
        #     continue
        # elif i > 1:
        #     break


        ra, dec, z = float(data[obj]['ra']), float(data[obj]['dec']), float(data[obj]['z'])

        print('\n[', list(data).index(obj)+1, '/', len(data), ']', obj, '|', data[obj]['ra'], data[obj]['dec'], data[obj]['z'])
        print('---------------------------------------------------------------------------------------------------------')

        transient_position = SkyCoord(ra, dec, unit=u.deg)
        try:
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                host_data = getTransientHosts(transientCoord=[transient_position], transientName=[obj], verbose=False, starcut="gentle", savepath=save_loc+'ghost_stuff/', GHOSTpath=GHOST_DATA)

            # If GLADE can't get magnitudes -- using NED name to get SDSS data
            print('Identified Host Galaxy:', host_data.loc[0, 'NED_name'])
            if np.isnan(host_data.loc[0, 'gKronMag']) or np.isnan(host_data.loc[0, 'iKronMag']):
                print('GLADE does not contain the g-mag & i-mag, searching SDSS...')
                host_position = SkyCoord(host_data['raMean'], host_data['decMean'], unit=u.deg)
                result = SDSS.query_crossid(host_position, photoobj_fields=['modelMag_g', 'modelMagErr_g', 'modelMag_i', 'modelMagErr_i'])

                # If GLADE or SDSS has magnitudes -- using PanSTARRS
                if result == None:
                    print('GLADE & SDSS do not contain the g-mag & i-mag, searching PanSTARRS...')
                    catalog_data = Catalogs.query_region(transient_position, radius=2 * u.arcsec, catalog="Panstarrs")
                    # print(catalog_data['rMeanKronMag'], catalog_data['iMeanKronMag'])

                    if len(catalog_data) == 0 or isinstance(catalog_data['gMeanKronMagErr'].value[0], np.ma.core.MaskedConstant) or isinstance(catalog_data['iMeanKronMagErr'].value[0], np.ma.core.MaskedConstant):
                        print('[!!!!!] GLADE, SDSS, and PanSTARRS do not contain the g-mag & i-mag data...')
                        gMag, iMag, iAbsMag = np.nan, np.nan, np.nan
                        gMagErr, iMagErr, iAbsMagErr = np.nan, np.nan, np.nan
                    else:
                        gMag, iMag, iAbsMag = catalog_data['gMeanKronMag'].value[0], catalog_data['gMeanKronMag'].value[0], (catalog_data['iMeanKronMag'].value[0] - cosmo.distmod(z).value)
                        gMagErr, iMagErr, iAbsMagErr = catalog_data['gMeanKronMagErr'].value[0], catalog_data['iMeanKronMagErr'].value[0], catalog_data['iMeanKronMagErr'].value[0]
                else:
                    # SDSS Results
                    gMag, iMag, iAbsMag = result['modelMag_g'].value[0], result['modelMag_i'].value[0], (result['modelMag_i'].value[0] - cosmo.distmod(z).value)
                    gMagErr, iMagErr, iAbsMagErr = result['modelMagErr_g'].value[0], result['modelMagErr_i'].value[0], result['modelMagErr_i'].value[0]
            else:
                # GLADE results
                gMag, iMag, iAbsMag = host_data['gKronMag'].loc[0], host_data['iKronMag'].loc[0], (host_data['iKronMag'].loc[0] - cosmo.distmod(z).value)
                gMagErr, iMagErr, iAbsMagErr = host_data['gKronMagErr'].loc[0], host_data['iKronMagErr'].loc[0], host_data['iKronMagErr'].loc[0]




            #  Mass Calculation -- Taylor et. al. 2011 -- eq. 8
            logstellarmass = (1.15 + (0.7*(gMag - iMag)) - (0.4*(iAbsMag)))

            # Error Propogation
            giMagErr = np.sqrt((gMagErr ** 2) + (iMagErr ** 2))
            logstellarmasserr = np.sqrt(((0.7**2)*(giMagErr**2)) + ((0.4**2)*(iAbsMagErr**2)))


        except Exception as error:
            print(obj+':', error)
            logstellarmass, logstellarmasserr = 0.00, 0.00

        if np.isnan(logstellarmass) == False and logstellarmass > 0.00:
            print('Success!', obj, 'host galaxy has a mass of:', logstellarmass, '+/-', logstellarmasserr, 'logM_* / [Msun]')
            all_z.append(z)
            all_logstellarmass.append(logstellarmass)
            if update_saved:
                data[obj].update({'logstellarmass': logstellarmass, 'logstellarmasserr': logstellarmasserr})
        else:
            print('Failed to find host galaxy!')
            if update_saved:
                data[obj].update({'logstellarmass': 0.00, 'logstellarmasserr': 0.00})

        # if i > 0:
        #     break
        # else:
        #     i =+ 1


    print('\nSuccessfully found mass of', len(all_z), '/', len(data), 'host galaxies!')
    if not keep_data:
        print('Removing GHOST data...')
        shutil.rmtree(save_loc+'ghost_stuff/') # Clear messy data
        os.mkdir(save_loc+'ghost_stuff/')
    if update_saved:
        print('Saving data to'+dict_path+'...')
        dict_handler(choice='pack', data_dict=data, path=dict_path)

    return all_z, all_logstellarmass
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

    failed_host_masses = []
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
    print('Failed host mass calculations:')
    print('[objname, ra, dec, z]')
    for fail in failed_host_masses:
        print(fail)

    if not keep_data:
        print('Removing GHOST data...')
        shutil.rmtree(save_loc + 'ghost_stuff/')  # Clear messy data
        os.mkdir(save_loc + 'ghost_stuff/')
    if update_saved:
        print('Saving data to' + dict_path + '...')
        dict_handler(choice='pack', data_dict=data, path=dict_path)
    return
# ==================================================================================================================== #
def salt3_fit(objs, plot_save_loc='../default/', plot_data=True, save_plot=True):
    print('[+++] Fitting data with SALT3...')

    CONSTANTS = get_constants()
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
            result, fitted_model = sncosmo.fit_lc(data, model, ['t0', 'x0', 'x1', 'c']) #  bounds={'z': (0.3, 0.7)}

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
            systime.sleep(3)
        except Exception as error:
            print(error)

    print('Successfully fit [', len(params), '/', len(objs), '] !')

    return params
# ==================================================================================================================== #
def lc_plot(objs, y_type = 'flux', pause_time=2, color_wheel = ['orange', 'cyan', 'violet', 'red', 'blue'],
            quiet=False, save_plots=True, save_loc='../snpy/misc_plots/'):
    print('[+++] Plotting LC data...')

    # Check quiet
    if quiet:
        sys.stdout = open(os.devnull, 'w')

    for obj in objs:
        print('-------------------------------------------------------------------------------------------------------')
        print('[', list(objs).index(obj)+1, '/', len(objs), ']')

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
        fig.suptitle("Hubble Residuals vs. " + x_params[1] + " of '"+(path.split('/')[-1][:-10]).upper()+"' 91bg-like SNe Ia\n" +  # Figure Title
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
def histogram_plotter(data_set='combined', raw=False, save_plot=False, save_loc='../default/', ignore_type=[], hist_bins = [40, 500, 40, 500, 80]):
    fig, ax = plt.subplots(3, 2, figsize=(12, 8), layout='constrained')

    # Get reviewed fits
    with open(get_constants()['reviewed_fits_txt'], 'r') as f:
        reviewed_good_fits, reviewed_okay_fits, reviewed_bad_fits = [], [], []
        for container in [reviewed_good_fits, reviewed_okay_fits, reviewed_bad_fits]:
            for i in f.readline().split(', ')[1:]:
                container.append(i[1:-1])
            container[-1] = container[-1][:-1]

    # Pull data
    data = np.genfromtxt(get_constants()[data_set+'_saved_loc']+data_set+'_saved.txt', dtype=str, skip_header=1, delimiter=', ')
    plot_titles = ['mu', 'st', 'Tmax', 'EBVhost', 'Host Mass']
    hist_indexs = [7, 8, 9, 10, 15]
    hist_err_indexs = [11, 12, 13, 14, 16]
    # xlimits = [None, [2.4e6, 2.5e6], None, [-0.5, 3], None]
    xlimits = [None, None, None, None, None]

    mu_hist, st_hist, Tmax_hist, EBV_hist, mass_hist = [], [], [], [], [] # 1, 7, 8, 9, 10, 16
    for i in range(len(data)):
        # Clean data
        if not raw:
            objname = data[i, 1]
            if objname[:2] == 'SN':
                objname = objname[2:]
            if ('good' in ignore_type) and (objname in reviewed_good_fits):
                continue
            if ('okay' in ignore_type) and (objname in reviewed_okay_fits):
                continue
            if ('bad' in ignore_type) and (objname in reviewed_bad_fits):
                continue

        mu_hist.append(float(data[i, 7]))
        st_hist.append(float(data[i, 8]))
        Tmax_hist.append(float(data[i, 9]))
        EBV_hist.append(float(data[i, 10]))
        mass_hist.append(float(data[i, 16]))


    # for arr in [mu_hist, st_hist, Tmax_hist, EBV_hist, mass_hist]:
    #     print(arr)

    arr_ind, arrs = 0, [mu_hist, st_hist, Tmax_hist, EBV_hist, mass_hist]
    for i in range(3):
        for j in range(2):
            ax[i, j].hist(arrs[arr_ind])
            arr_ind += 1
            print(arr_ind)
            # plt.suptitle("Parameters for '" + data_set.upper() + "' data\n Number of Transients: " + str(len(objs[:, 0])), fontsize=20)
            # ax[i, j].hist(objs[:, hist_indexs[order]].astype(float), bins=hist_bins[order])
            # ax[i, j].set_title(plot_titles[order])
            # ax[i, j].set_xlim(xlimits[order])
    plt.show()




    # order = 0
    # for j in range(size[1]):
    #     i = 0
    #     while i < size[0]:
    #         plt.suptitle("Parameters for '"+data_set.upper()+"' data\n Number of Transients: " + str(len(objs[:, 0])), fontsize=20)
    #         ax[i, j].hist(objs[:, hist_indexs[order]].astype(float), bins=hist_bins[order])
    #         ax[i, j].set_title(plot_titles[order])
    #         ax[i, j].set_xlim(xlimits[order])
    #
    #         i += 1
    #         order += 1
    #         if order > len(plot_titles)-2: # 2 because I havent ran ghost
    #             break
    # if save_plot:
    #     plt.savefig(save_loc+'hist_'+data_set.upper()+'.png')
    # plt.show()

    # # Cleaning stretch data
    # true_vals = np.where(objs[:, 8].astype(float)-2400000 > 0)[0]
    # print(objs[:, 8].astype(float)[true_vals])
    # plt.hist(objs[:, 8].astype(float)[true_vals])
    # plt.show()

    return