import os
import sys
import snpy
import glob
import time as systime
import numpy as np
import matplotlib.pyplot as plt
from zipfile import ZipFile
from scripts import tns_redshifts

import scripts.general as gen

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
def TNS_unpack(obj_ra, obj_dec):
    print('Searching TNS for ['+str(obj_ra)+', '+str(obj_dec)+']...')
    tns_key = np.genfromtxt('../working_data/TNS_key_new.txt', dtype=str, delimiter=', ', skip_header=1)

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
        details = gen.TNS_details(obj_ra, obj_dec)
        obj_name, obj_z = details['objname'], details['redshift']
        print('Found! ['+str(obj_ra)+', '+str(obj_dec)+'] -- ', obj_name, ',', obj_z)
        # Save to key
        with open('../working_data/TNS_key_new.txt', 'a') as f:
            f.write(str(obj_ra) + ', ' + str(obj_dec) + ', ' + obj_name + ', ' + str(obj_z) + '\n')

    return obj_name, obj_z
def newer_snpy_fit(path, objname, save_loc, plot_save_loc,
                 skip_problems=False, use_saved=False, snpy_plots=True, save_plots=True):
    problem_children = handle_problem_children(state='READ') # Open problem children

    # Load data
    n_s = snpy.get_sn(path)
    n_s.choose_model('EBV_model2', stype='st')
    n_s.set_restbands()  # Auto pick appropriate rest-bands

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
            n_s.save(save_loc + objname + '_EBV_model2.snpy')
            run = False

            if snpy_plots:
                n_s.plot(outfile=plot_save_loc + objname + '_snpyplots.png')
                plt.show()
            plt.close()

            return {'ra': n_s.ra, 'dec': n_s.decl, 'z': n_s.z, 'MJDs': mjds, 'MJDe': mjde,
                    'mu': n_s.parameters['DM'], 'st': n_s.parameters['st'], 'Tmax': n_s.parameters['Tmax'],
                    'EBVhost': n_s.parameters['EBVhost'],
                    'mu_err': n_s.errors['DM'], 'st_err': n_s.errors['st'], 'Tmax_err': n_s.errors['Tmax'],
                    'EBVhost_err': n_s.errors['EBVhost']}

        except Exception as error:
            if 'All weights for filter' and 'are zero.' in str(error):
                print('[!!!] Weights for filter', str(error).split(' ')[4], 'are zero. Removing...')
                del n_s.data[str(error).split(' ')[4]]
            elif 'Error:  to solve for EBVhost, you need to fit more than one filter' in str(error):
                # print('Only one filter, cannot fit data.')
                problem_children = np.append(problem_children, objname)
                handle_problem_children(state='WRITE', problem_c=problem_children) # Commit problem children
                return error
            else:
                problem_children = np.append(problem_children, objname)
                handle_problem_children(state='WRITE', problem_c=problem_children) # Commit problem children
                return error
def snpy_fit_indv(objname):
    path = '/content/HiloCATsSN1991bg/snpy/atlas/'+objname+'_EBV_model2.snpy'

    n_s = snpy.get_sn(path)
    n_s.choose_model('EBV_model2', stype='st')
    n_s.set_restbands() # Auto pick appropriate rest-bands

    # Sort out empty filters & get start and end time
    mjds, mjde = [], []
    filter_wheel = []
    for filter in list(n_s.data.keys()):
        if len(n_s.data[filter].MJD) <= 3:
            print('\t', objname, 'has too few points in', filter, 'filter')
            continue
        mjds.append(min(n_s.data[filter].MJD))
        mjde.append(max(n_s.data[filter].MJD))
        filter_wheel.append(filter)

    n_s.fit(bands=filter_wheel, dokcorr=True, k_stretch=False, reset_kcorrs=True, **{'mangle':1,'calibration':0})
    plt.show()

    print('Results:',
        'mu =', n_s.parameters['DM'], '+/-', n_s.errors['DM'],'\n',
        '\t st =', n_s.parameters['st'], '+/-', n_s.errors['st'],'\n',
        '\t Tmax =', n_s.parameters['Tmax'], '+/-', n_s.errors['Tmax'],'\n',
        '\t EBVhost =', n_s.parameters['EBVhost'],  '+/-', n_s.errors['EBVhost'],'\n',
        '\t MJD min:', min(mjds), '| MJD max:', max(mjde))

    return
def read_DR3(loc='../txts/DR3_fits.dat'):
    data = np.genfromtxt(loc, dtype=str, skip_header=1)
    dr3 = {}
    for n in range(len(data[:, 0])):
        dr3.update({data[:, 0][n]: {'st': float(data[:, 1][n]), 'e_st': float(data[:, 2][n]), 'z': float(data[:, 3][n]),
                           'Tmax': float(data[:, 5][n]), 'e_Tmax': float(data[:, 6][n]),
                           'EBVHost': float(data[:, 25][n]), 'e_EBVHost': float(data[:, 26][n])}})
    return dr3
def dict_unpacker(path, delimiter=', '):
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
def dict_packer(data_dict, save_loc, delimiter=', '):
    print('[+++] Packing objects into text file...')

    catagories = list(data_dict[list(data_dict.keys())[0]].keys())
    with open(save_loc, 'w') as f:
        f.write('objname')
        for category in catagories:
            f.write(delimiter+category)
        f.write('\n')
        for objname in data_dict:
            f.write(objname)
            for category in catagories:
                f.write(delimiter+str(data_dict[objname][category]))
            f.write('\n')
    print('[+++] Files packed to', save_loc)
    return
def handle_problem_children(state, problem_c=None):
    path = '../problem_children.txt'

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
def save_to_zip(zip_loc, save_loc):
    print('Saving zipped files to...', save_loc)
    files = glob.glob(zip_loc+'*.png')
    with ZipFile(save_loc, 'w') as zip:
        for n_file in files:
            zip.write(n_file)
    return
# =================================================================================================================== #
def data_proccesser(data_set = 'ZTF', mag_unc_max=1, quiet=False):
    print('[+++] Processing '+data_set+' data...')

    # Check quiet
    if quiet:
        sys.stdout = open(os.devnull, 'w')

    # Chose data set
    if data_set == 'ATLAS':
        hdr_len = 1
        data_delim, hdr_delim = ',', ','
        files = glob.glob('../data/ATLAS/*.txt')
        d_i = {'mag': 3, 'dmag': 4, 'filter': 6, 'time': 8, 'flux': 16, 'dflux': 17, 'uJy': 24, 'duJy': 25}
    elif data_set == 'ZTF':
        hdr_len = 56
        data_delim, hdr_delim = None, ', '
        files = glob.glob('../data/ZTF/*.txt')
        d_i = {'filter': 4, 'time': 22, 'flux': 24, 'dflux': 25}
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

        # Object Creation
        obj_filters, obj_time, obj_flux, obj_flux_unc, obj_mag, obj_mag_unc = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
        for i in range(len(data[:, 0])):
            # Flux to magnitude
            if data_set == 'ZTF':
                if data[i, d_i['flux']] == 'null':
                    continue
                else:
                    n_flux, n_flux_unc = float(data[i, d_i['flux']]), float(data[i, d_i['dflux']])
                    n_mag = -2.5 * np.log10(n_flux) + 23.9
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
            if n_flux <= 0: # Negatives
                continue
            if n_mag <= 0:  # Negatives
                continue
            if n_mag_unc <= 0:  # Negatives
                continue
            if n_mag_unc > mag_unc_max:
                continue

            # Commit variables
            obj_filters = np.append(obj_filters, data[i, d_i['filter']])
            obj_time = np.append(obj_time, float(data[i, d_i['time']]))
            obj_flux = np.append(obj_flux, n_flux)
            obj_flux_unc = np.append(obj_flux_unc, n_flux_unc)
            obj_mag = np.append(obj_mag, n_mag)
            obj_mag_unc = np.append(obj_mag_unc, n_mag_unc)

        # Create Object
        obj = {'ra': obj_ra, 'dec': obj_dec, 'z': obj_z, 'filters': obj_filters, 'time': obj_time,
               'flux': obj_flux, 'dflux': obj_flux_unc, 'mag': obj_mag, 'dmag': obj_mag_unc}
        print('Retrived: ' + obj_name + ' | ra: '+str(obj_ra)+' \ dec: '+str(obj_dec))

        # Commit Object
        data_objs.update({obj_name: obj})

    # Restore print statements
    sys.stdout = sys.__stdout__

    return data_objs
def write_ASCII(objs, data_set, save_loc):
    print('[+++] Saving data to ASCII files for SNooPy...')
    for obj in objs:
        with open(save_loc + obj + '_snpy.txt', 'w') as f:
            # Get list of all filters
            if data_set == 'ATLAS':
                filter_wheel = ['c', 'o']
                filter_label = ['ATgr', 'ATri']
            elif data_set == 'ZTF':
                filter_wheel = ['ZTF_g', 'ZTF_r', 'ZTF_i']
                filter_label = ['g', 'r', 'i']

            # File Header
            # Line 1 -- Objname, Helio-Z, RA, Dec (Ex. SN1981D 0.005871 50.65992 -37.23272)
            f.write(str(obj)+' '+str(objs[obj]['z'])+' '+str(objs[obj]['ra'])+' '+str(objs[obj]['dec'])+'\n')

            for f_w in filter_wheel:
                f_indexs = np.where(objs[obj]['filters'] == f_w)[0]

                # 'o'/'ATri'-filter photometry block -- Date (JD/MJD), mag, err (674.8593 12.94 0.11)
                # 'c'/'ATgr'-filter photometry block
                f.write('filter '+filter_label[filter_wheel.index(f_w)]+'\n')
                for i in range(len(objs[obj]['time'][f_indexs])):
                    f.write(str(objs[obj]['time'][i])+'\t'+str(objs[obj]['mag'][i])+'\t'+str(objs[obj]['dmag'][i])+'\n')
    return
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
            for saved_fit in glob.glob(save_loc+'*.snpy'):
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
        n_s = snpy.get_sn(path)
        n_s.choose_model('EBV_model2', stype='st')
        n_s.set_restbands()  # Auto pick appropriate rest-bands

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
                n_s.save(save_loc + objname + '_EBV_model2.snpy')
                run = False

                if snpy_plots:
                    n_s.plot(outfile=save_loc+'plots/' + objname + '_snpyplots.png')
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
