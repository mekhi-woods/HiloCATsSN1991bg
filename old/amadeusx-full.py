""" IMPORTS """
import os
import glob
import json
import snpy
import scipy
import shutil
import sncosmo
import requests
import urllib.request
import time as systime
import numpy as np
import matplotlib.pyplot as plt
from zipfile import ZipFile
from requests.auth import HTTPBasicAuth
from HiloCATsSN1991bg.scripts import tns_redshifts
import astro_ghost
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from astropy import cosmology as cosmo
from astro_ghost.ghostHelperFunctions import getTransientHosts

""" GLOBALS """
# TNS CREDINTIALS
tns_bot_id, tns_bot_name, tns_bot_api_key = '73181', 'YSE_Bot1', '0d771345fa6b876a5bb99cd5042ab8b5ae91fc67'

# PATHS
ROOT_PATH = '/content/HiloCATsSN1991bg/'
SNPY_ROOT = ROOT_PATH+'snpy/'
PLOTS_ROOT = ROOT_PATH+'plots/'
TEST_ROOT = ROOT_PATH+'test/'

SNPY_BURNS = SNPY_ROOT+'burns/'
SNPY_BURNS_PLOTS = SNPY_BURNS+'plots/'
BURNS_SAVE_TXT = SNPY_BURNS+'burns_saved.txt'

DATA_ATLAS = ROOT_PATH+'data/ATLAS/'
SNPY_ATLAS = SNPY_ROOT+'atlas/'
SNPY_ATLAS_PLOTS = SNPY_ATLAS+'plots/'
SNPY_ATLAS_ASCII = SNPY_ATLAS+'ascii/'
ATLAS_SAVE_TXT = SNPY_ATLAS+'atlas_saved.txt'

DATA_ZTF = ROOT_PATH+'data/ZTF/'

GHOST_DATA = ROOT_PATH+'data/GHOST/'

PROB_CHILD_TXT = ROOT_PATH+'problem_children.txt'
TNS_KEY_TXT = ROOT_PATH+'TNS_key.txt'

COSMO_MODEL = cosmo.FlatLambdaCDM(H0=70, Om0=0.3)

""" ZTF """
def ztf_collection(submit=False, limit=1000):
    objs = dict_unpacker(ATLAS_SAVE_TXT)
    print("Number of targets =", len(objs))
    ralist, declist, jdslist, jdelist = [], [], [], []
    i = 0
    for obj in objs:
        ralist.append(float(objs[obj]['ra']))
        declist.append(float(objs[obj]['dec']))
        jdslist.append(float(objs[obj]['MJDs'])-100)
        jdelist.append(float(objs[obj]['MJDe'])+100)
        i += 1
        if i % limit == 0:
            if submit:
                ztf_submit_post(ralist, declist, jdslist, jdelist)
            ralist, declist, jdslist, jdelist = [], [], [], []
    if submit and len(ralist) > 0:
        ztf_submit_post(ralist, declist, jdslist, jdelist)

    return

def ztf_submit_post(ra_list, dec_list, jds, jde):
    print('Submiting request to ZTF...')

    email = 'mekhidw@hawaii.edu' # email you subscribed with.
    userpass = 'wxdk286' # password that was issued to you.

    ra, dec = json.dumps(ra_list), json.dumps(dec_list)
    jdend, jdstart = json.dumps(jde), json.dumps(jds)
    payload = {'ra': ra, 'dec': dec, 'jdstart': jdstart, 'jdend': jdend, 'email': email, 'userpass': userpass}

    # fixed IP address/URL where requests are submitted:
    url = 'https://ztfweb.ipac.caltech.edu/cgi-bin/batchfp.py/submit'
    r = requests.post(url, auth=('ztffps', 'dontgocrazy!'), data=payload)

    print('RA ['+str(type(payload['ra'][0]))+']:', payload['ra'])
    print('DEC ['+str(type(payload['dec'][0]))+']:', payload['dec'])
    print('MJD Start ['+str(type(payload['jdstart'][0]))+']:', payload['jdstart'])
    print('MJD End ['+str(type(payload['jdend'][0]))+']:', payload['jdend'])
    print("Status_code=", r.status_code)
    return

def ztf_alt_collection():
    ra, dec, jds, jde = 186.07860833333334, 10.446222222222222, 2460350, 2460450

    email = 'mekhidw@hawaii.edu' # email you subscribed with.
    userpass = 'wxdk286' # password that was issued to you.

    cmd = f"wget --http-user=ztffps --http-passwd=dontgocrazy! -O log.txt \"https://ztfweb.ipac.caltech.edu/cgi-bin/requestForcedPhotometry.cgi?ra={dec}&dec={ra}&jdstart={jds}&jdend={jde}&email={email}&userpass={userpass}\""
    print(cmd)
    os.system(cmd)

    return

""" ATLAS """
def atlas_collection(quiet=False, check_data=True):
    api_token = '7f4e1dee8f52cf0c8cdaf55bf29d04bef4810fb4'

    if check_data and len(glob.glob(DATA_ATLAS+'*.txt')) > 1:
        print('[+++] ATLAS data already collected, passing step...')
        return
    else:
        print('[+++] No data detected, collecting ATLAS data...')
        if os.path.exists(DATA_ATLAS+'tmp.npz'):
            pickle = np.load(DATA_ATLAS+'tmp.npz', allow_pickle=True)
            data = pickle['data']

        data = requests.post('https://star.pst.qub.ac.uk/sne/atlas4/api/objectlist/',
                             headers={'Authorization': f'Token {api_token}'},
                             data={'objectlistid':2}).json()

        np.savez(DATA_ATLAS+'tmp.npz', data=data)

        count = 0
        for d in data:
            if d['observation_status'] is not None and d['observation_status'].startswith('SN Ia') and '91bg' in d['observation_status']:
                count += 1
                if not quiet:
                    print(d['atlas_designation'],d['observation_status'].replace(' ',''),d['ra'],d['dec'])


                ids = d['id']
                base_url = 'https://star.pst.qub.ac.uk/sne/atlas4/lightcurveforced/1161048951013729300/'
                new_url = base_url.replace('1161048951013729300/',str(ids))
                if not quiet:
                    print(new_url)

                idfile = DATA_ATLAS+'/' + str(ids)+'.txt'
                if os.path.exists(idfile):
                    continue
                urllib.request.urlretrieve(str(new_url), str(idfile))
                if not quiet:
                    print(idfile)

            if count > 300:
                break
    return

def atlas_tns_collection(quiet=False, t_sleep=5, max_attempts=5):
    files = glob.glob(DATA_ATLAS+'/*.txt')
    tns_key = dict_unpacker(TNS_KEY_TXT)

    attempts = max_attempts
    while attempts > 0:
        try:
            f = open(TNS_KEY_TXT, 'a')
            if len(tns_key) == 0:
                f.write('atlas_name, objname, ra, dec, z\n')
            for n in range(len(files)):
                tracker = '['+str(n+1)+'/'+str(len(files))+']' # Purely cosmetic
                atlas_name = files[n][len(DATA_ATLAS):-4]
                data = np.genfromtxt(files[n], dtype=float, delimiter=',', skip_header=1)
                ra = np.average(data[:, 1])
                dec = np.average(data[:, 2])

                if atlas_name not in tns_key:
                    details = TNS_details(ra, dec)
                    tns_key.update({atlas_name: {'objname': details['objname'], 'ra': ra, 'dec': dec, 'z': details['redshift']}})
                    f.write(atlas_name+', '+str(details['objname'])+', '+str(ra)+', '+str(dec)+', '+str(details['redshift'])+'\n')
                    if not quiet:
                        print(tracker, details['objname']+':', ra, '|', dec, '|', details['redshift'])
            f.close()
        except:
            print('[', attempts, '/', max_attempts, '] Hit time error, sleeping for '+str(t_sleep)+' seconds...')
            systime.sleep(t_sleep)
            attempts -= 1
    if attempts == 0:
        print('[!!!] ATTEMPT LIMIT HIT')
        return

    dict_packer(tns_key, TNS_KEY_TXT)

    return tns_key

def atlas_processing(err_max=100, n_iter=10, sleep_t=5, use_TNS=False):
    print('[+++] Retrieving data from...', DATA_ATLAS)
    # Retrieve file paths
    files = glob.glob(DATA_ATLAS+'/*.txt')
    if n_iter != 0 and n_iter <= len(files):
        files = files[:n_iter]

    objs = {}
    for n in range(len(files)):
        # Tracking/Cosmetics
        ATLAS_name = files[n][len(DATA_ATLAS):-4]
        tracker = '['+str(n+1)+'/'+str(len(files))+']' # Purely cosmetic
        print(tracker, '\n', '\t\tPulling', ATLAS_name, 'data...')

        # Reading file path
        try:
            data = np.genfromtxt(files[n], dtype=str, delimiter=',', skip_header=1)
            if len(data) == 0: # Handling empty files
                print('[!!!] \t\tFile '+ATLAS_name+' empty...skipping') # If empty, skips
                continue
        except:
            print('[!!!] \t\tUnknown error, skipping...') # Numpy was doing a weird thing, so crash hander until i figure it out
            continue
        ra, dec = np.average(data[:, 1].astype(float)), np.average(data[:, 2].astype(float)) # Recoring RA & DEC (average of columns)

        # Using TNS to get object name
        tns_key = dict_unpacker(TNS_KEY_TXT)
        if use_TNS and (ATLAS_name in tns_key):
            print('\t\tAttempting to use TNS key for objname...')
            objname = tns_key[ATLAS_name]['objname']
            z = tns_key[ATLAS_name]['z']
        else:
            print('\t\tRetrieving TNS data for...', ATLAS_name, '[sleeping for', sleep_t, 'seconds...]')
            systime.sleep(sleep_t)
            details = TNS_details(ra, dec)
            objname = details['objname']
            z = details['redshift']
            with open(TNS_KEY_TXT, 'a') as f:
                f.write(ATLAS_name+', '+str(objname)+', '+str(ra)+', '+str(dec)+', '+str(z)+'\n')
            print('\t\tSaved TNS data to TNS key...')

        mag = np.char.replace(data[:, 3], '>', '') # Removes greater than symbols
        dmag, filters, time, flux, dflux = data[:, 4], data[:, 6], data[:, 8], data[:, 24], data[:, 25] # Reads rest of categories
        objs.update({objname: {'ra': ra, 'dec': dec, 'z': z, 'time': time, 'flux': flux, 'dflux': dflux, 'mag': mag, 'dmag': dmag, 'filters': filters}})

        ## SLICING DATA
        # ------------------------------------------------------------------------------------------------------------------------------------------------------
        # Remove 'None' from categories
        mod_empty = np.array([])
        for cat in ['time', 'flux', 'dflux', 'mag', 'dmag', 'filters']:
            mod_empty = np.append(mod_empty, np.where(objs[objname][cat] == 'None')[0])
        mod_empty = np.unique(mod_empty).astype(int) # Remove dupes
        for cat in ['time', 'flux', 'dflux', 'mag', 'dmag', 'filters']:
            objs[objname][cat] = np.delete(objs[objname][cat], mod_empty)

        # Finds negative fluxes
        mod_positive = np.array([])
        for cat in ['time', 'flux', 'dflux', 'mag', 'dmag']:
            mod_positive = np.append(mod_positive, np.where(objs[objname][cat].astype(float) <= 0)[0])
        mod_positive = np.unique(mod_positive).astype(int) # Remove dupes
        for cat in ['time', 'flux', 'dflux', 'mag', 'dmag', 'filters']:
            objs[objname][cat] = np.delete(objs[objname][cat], mod_positive)

        # Find outliers beyond error limit
        mod_err = np.array([])
        for cat in ['dflux', 'dmag']:
            mod_err = np.append(mod_err, np.where(np.abs(objs[objname][cat].astype(float)) > err_max)[0])
        mod_err = np.unique(mod_err).astype(int) # Remove dupes
        for cat in ['time', 'flux', 'dflux', 'mag', 'dmag', 'filters']:
            objs[objname][cat] = np.delete(objs[objname][cat], mod_err)

        # Set everything as floats
        for cat in ['time', 'flux', 'dflux', 'mag', 'dmag']:
            objs[objname][cat] = objs[objname][cat].astype(float)

        # Seperate into orange & cyan channels
        for cat in ['time', 'flux', 'dflux', 'mag', 'dmag']:
            objs[objname].update({cat+'_o': objs[objname][cat][np.where(objs[objname]['filters'] == 'o')[0]]})
            objs[objname].update({cat+'_c': objs[objname][cat][np.where(objs[objname]['filters'] == 'c')[0]]})
        # ------------------------------------------------------------------------------------------------------------------------------------------------------

        print('\t\tRetrived:', objname+' |', 'ra:', ra, '\\', 'dec:', dec)
    print('[!!!]\t\tRetrieved & processed', len(objs), 'SNe from ATLAS!')
    return objs

def atlas_write_ASCII(atlas_objs, save_loc, quiet=True):
    print('[+++] Saving data to ASCII files for SNooPy...')
    for obj in atlas_objs:
        with open(save_loc+obj+'_snpy.txt', 'w') as f:
            # Line 1 -- Objname, Helio-Z, RA, Dec (Ex. SN1981D 0.005871 50.65992 -37.23272)
            f.write(str(obj)+' '+str(atlas_objs[obj]['z'])+' '+str(atlas_objs[obj]['ra'])+' '+str(atlas_objs[obj]['dec'])+'\n')

            # 'o'/'ATri'-filter photometry block -- Date (JD/MJD), mag, err (674.8593 12.94 0.11)
            f.write('filter ATri\n')
            for i in range(len(atlas_objs[obj]['time_o'])):
                f.write(str(atlas_objs[obj]['time_o'][i])+'\t'+str(atlas_objs[obj]['mag_o'][i])+'\t'+str(atlas_objs[obj]['dmag_o'][i])+'\n')

            # # 'c'/'ATgr'-filter photometry block
            f.write('filter ATgr\n')
            for i in range(len(atlas_objs[obj]['time_c'])):
                f.write(str(atlas_objs[obj]['time_c'][i])+'\t'+str(atlas_objs[obj]['mag_c'][i])+'\t'+str(atlas_objs[obj]['dmag_c'][i])+'\n')
    return

def atlas_snpy_fitting(n_iter=0, skip_problems=True, use_saved=True, snpy_plots=True, save_plots=True):
    print('[+++] Fitting ATLAS data with SNooPy...')
    fit_args = {'skip_problems': skip_problems, 'use_saved': use_saved, 'snpy_plots': snpy_plots, 'save_plots': save_plots}
    print('Fitting arguments: ', fit_args)
    objpaths = glob.glob(SNPY_ATLAS_ASCII+'*')
    if n_iter != 0 and n_iter <= len(objpaths):
        objpaths = objpaths[:n_iter]

    objs = {}
    err_i = 0
    for n in range(len(objpaths)):
        tracker = '['+str(n+1)+'/'+str(len(objpaths))+']' # Purely cosmetic
        objname = objpaths[n][len(SNPY_ATLAS_ASCII):-9]
        print(tracker, objname)
        temp_dict = snpy_fit(objpaths[n], objname, save_loc=SNPY_ATLAS, plot_save_loc=SNPY_ATLAS_PLOTS, **fit_args)

        if type(temp_dict) == dict:
            print('\tResults: \n',
                  '\t\tmu =', temp_dict['mu'], '+/-', temp_dict['mu_err'],'\n',
                  '\t\tst =', temp_dict['st'], '+/-', temp_dict['st_err'],'\n',
                  '\t\tTmax =', temp_dict['Tmax'], '+/-', temp_dict['Tmax_err'],'\n',
                  '\t\tEBVhost =', temp_dict['EBVhost'],  '+/-', temp_dict['EBVhost_err'])
            print('\tMJD min:', temp_dict['MJDs'], '+/- 100 MJD', '| MJD max:', temp_dict['MJDe'], '+/- 100 MJD\n')
            objs.update({objname: temp_dict})
        else:
            err_i += 1
            print('[!!!]\t', temp_dict, '\n')
    print('Finshed! Successfully fit', len(objpaths)-err_i, '/', len(objpaths), 'SNe from ATLAS! ['+str(round(((len(objpaths)-err_i) / len(objpaths))*100, 2))+'%]')

    # Save Data
    if len(objs) > 0:
        dict_packer(objs, ATLAS_SAVE_TXT, delimiter=', ') # Save data from fitting
    if snpy_plots and save_plots:
        save_to_zip(SNPY_ATLAS_PLOTS, SNPY_ATLAS_PLOTS+'atlas_snpy_plots.zip')

    return

def atlas_plotting(choice):
    if 'reg_hist' in choice:
        print('[+++] Ploting Histogram of ATLAS SNooPy Fitting Parameters...')
        data = np.genfromtxt(SNPY_ATLAS+'atlas_saved.txt', dtype=str, delimiter=', ', skip_header=1)
        objnames, st, Tmax, EBVHost = data[:, 0], data[:, 3].astype(float), data[:, 4].astype(float), data[:, 5].astype(float)
        reg_params = [st, Tmax, EBVHost]
        bins_reg = [15, 20, 15]
        plot_title = 'ATLAS SNooPy Parameters'

        fig, ax = plt.subplots(1, 3, figsize=(15, 5))
        titles = ['st', plot_title+'\nTmax', 'EBVhost']
        for i in range(len(titles)):
            ax[i].hist(reg_params[i], bins=bins_reg[i])
            ax[i].set_title(titles[i])
        plt.show()

""" BURNS/CSP """
def burns_cspvdr3(fit_filters=None, skip_problems=False, use_saved=False, snpy_plots=True):
    print('[+++] Fitting CSP data with SNooPy...')
    # Get Chris Burns Data
    objnames = np.genfromtxt('/content/HiloCATsSN1991bg/txts/91bglike_justnames.txt', dtype=str, delimiter=', ')

    # Get CSP paths of Chris Burns Data
    objpaths = []
    for name in objnames:
        if os.path.exists('/content/HiloCATsSN1991bg/data/CSPdata/SN'+name+'_snpy.txt'):
            objpaths.append('/content/HiloCATsSN1991bg/data/CSPdata/SN'+name+'_snpy.txt')
        else:
            print(name, 'not found...')

    # Fitting
    objs = {}
    err_i = 0
    fit_args = {'skip_problems': skip_problems, 'use_saved': use_saved, 'snpy_plots': snpy_plots}
    print('Fitting arguments: ', fit_args)
    for n in range(len(objpaths)):
        tracker = '['+str(n+1)+'/'+str(len(objpaths))+']' # Purely cosmetic
        objname = objpaths[n][39:-9]
        print(tracker, objname)
        temp_dict = snpy_fit(objpaths[n], objname, save_loc=SNPY_BURNS, plot_save_loc=SNPY_BURNS_PLOTS, **fit_args)

        if type(temp_dict) == dict:
            print('\tResults: \n',
                  '\t\tmu =', temp_dict['mu'], '+/-', temp_dict['mu_err'],'\n',
                  '\t\tst =', temp_dict['st'], '+/-', temp_dict['st_err'],'\n',
                  '\t\tTmax =', temp_dict['Tmax'], '+/-', temp_dict['Tmax_err'],'\n',
                  '\t\tEBVhost =', temp_dict['EBVhost'],  '+/-', temp_dict['EBVhost_err'])
            print('\tMJD min:', temp_dict['MJDs'], '+/- 100 MJD', '| MJD max:', temp_dict['MJDe'], '+/- 100 MJD\n')
            objs.update({objname: temp_dict})
        else:
            err_i += 1
            print('[!!!]\t', temp_dict, '\n')
    print('Finshed! Successfully fit', len(objpaths)-err_i, '/', len(objpaths), 'SNe from ATLAS! ['+str(round(((len(objpaths)-err_i) / len(objpaths))*100, 2))+'%]')

    # Announce problem children
    print('Problem children: ', handle_problem_children(state='READ'))

    # Save Data
    if len(objs) > 0:
        dict_packer(objs, BURNS_SAVE_TXT, delimiter=', ') # Save data from fitting
    return

def burns_plotting(choice):
    if 'reg_hist' in choice:
        print('[+++] Ploting Histogram of SNooPy Fitting Parameters...')
        data = np.genfromtxt(BURNS_SAVE_TXT, dtype=str, delimiter=', ', skip_header=1)
        objnames, st, Tmax, EBVHost = data[:, 0], data[:, 3].astype(float), data[:, 4].astype(float), data[:, 5].astype(float)
        reg_params = [st, Tmax, EBVHost]
        bins_reg = [15, 20, 15]
        plot_title = 'CSP SNooPy Parameters'

        fig, ax = plt.subplots(1, 3, figsize=(15, 5))
        titles = ['st', plot_title+'\nTmax', 'EBVhost']
        for i in range(len(titles)):
            ax[i].hist(reg_params[i], bins=bins_reg[i])
            ax[i].set_title(titles[i])
        plt.show()

    if 'res_hist' in choice:
        print('[+++] Ploting Histogram of SNooPy-Chris Burns Parameters Residuals...')
        data = np.genfromtxt(BURNS_SAVE_TXT, dtype=str, delimiter=', ', skip_header=1)
        objnames, st, Tmax, EBVHost = data[:, 0], data[:, 3].astype(float), data[:, 4].astype(float), data[:, 5].astype(float)

        # Take the difference between CSP and DR3 file
        st_res, Tmax_res, EBVHost_res = [], [], []
        dr3 = read_DR3()
        for n in range(len(objnames)):
            if objnames[n] in dr3:
                st_res.append(st[n] - dr3[objnames[n]]['st'])
                Tmax_res.append(Tmax[n] - dr3[objnames[n]]['Tmax'])
                EBVHost_res.append(EBVHost[n] - dr3[objnames[n]]['EBVHost'])\

        # Correct for MJD -- 53000
        for i in range(len(Tmax_res)):
            Tmax_res[i] = Tmax_res[i] + 53000

        # Plot
        res_params = [st_res, Tmax_res, EBVHost_res]
        bins_res = [30, 50, 20]
        xlims = [[-0.1, 0.1], [-3, 3], [-0.15, 0.15]]
        plot_title = 'SNooPy-Chris Burns Parameters'
        fig, ax = plt.subplots(1, 3, figsize=(15, 5))
        titles = ['st', plot_title+'\nTmax', 'EBVhost']
        for i in range(len(titles)):
            ax[i].hist(res_params[i], bins=bins_res[i])
            ax[i].set_title(titles[i])
            ax[i].set_xlim(xlims[i])

        plt.show()

    if 'zvmu' in choice:
        print('[+++] Ploting redshift [z] vs distance mod [mu] of SNooPy Parameters...')
        data = np.genfromtxt(BURNS_SAVE_TXT, dtype=str, delimiter=', ', skip_header=1)
        objnames, mu, z = data[:, 0], data[:, 1].astype(float), data[:, 2].astype(float)

        plt.figure(figsize=(8, 4))
        for n in range(len(objnames)):
            plt.loglog(z[n], mu[n], label=objnames[n], marker='o')
        plt.title('CSP Redshift vs. Distance Mod\n SNe # = '+str(len(objnames)))
        plt.xlabel('Redshift')
        plt.ylabel('Distance Mod')
        plt.legend()
        plt.show()

    return

""" GHOST """
def ghost_host_galaxy(dict_path, save_loc=TEST_ROOT, keep_data=True, update_saved=False):
    cosmo = FlatLambdaCDM(70, 0.3) # Hubble Constant, Omega-Matter
    data = dict_unpacker(dict_path)
    all_z, all_logstellarmass = [], []
    print('[+++] Finding host galaxy mass using GHOST...')

    i = 0
    for obj in data:
        ra, dec, z = float(data[obj]['ra']), float(data[obj]['dec']), float(data[obj]['z'])

        print('\n[', list(data).index(obj)+1, '/', len(data), ']', obj, '|', data[obj]['ra'], data[obj]['dec'], data[obj]['z'])
        print('---------------------------------------------------------------------------------------------------------')

        transient_position = SkyCoord(ra, dec, unit=u.deg)
        try:
            host_data = getTransientHosts(transientCoord=[transient_position], transientName=[obj], verbose=False, starcut="gentle", savepath=save_loc+'ghost_stuff/', GHOSTpath=GHOST_DATA)

            gMag, iMag, iAbsMag = host_data['gKronMag'].loc[0], host_data['iKronMag'].loc[0], (host_data['iKronMag'].loc[0] - cosmo.distmod(z).value)
            gMagErr, iMagErr, iAbsMagErr = host_data['gKronMagErr'].loc[0], host_data['iKronMagErr'].loc[0], host_data['iKronMagErr'].loc[0]
            giMagErr = np.sqrt((gMagErr**2) + (iMagErr**2))

            logstellarmass = (1.15 + (0.7*(gMag - iMag)) - (0.4*(iAbsMag))) # Taylor et. al. 2011 -- eq. 8
            logstellarmasserr = np.sqrt(((0.7**2)*(giMagErr**2)) + ((0.4**2)*(iAbsMagErr**2)))

        except Exception as error:
            print(obj+':', error)
            logstellarmass, logstellarmasserr = 0.00, 0.00

        if logstellarmass != np.nan and logstellarmass > 0:
            print('Success!', obj, 'host galaxy has a mass of:', logstellarmass, '+/-', logstellarmasserr, 'logM_* / [Msun]')
            all_z.append(z)
            all_logstellarmass.append(logstellarmass)
            if update_saved:
                data[obj].update({'logstellarmass': logstellarmass, 'logstellarmasserr': logstellarmasserr})
        else:
            print('Failed to find host galaxy!')
            if update_saved:
                data[obj].update({'logstellarmass': 0.00, 'logstellarmasserr': 0.00})

    print('\nSuccessfully found mass of', len(all_z), '/', len(data), 'host galaxies!')
    if not keep_data:
        print('Removing GHOST data...')
        shutil.rmtree(save_loc+'ghost_stuff/') # Clear messy data
    if update_saved:
        print('Saving data to'+dict_path+'...')
        dict_packer(data, dict_path)

    return all_z, all_logstellarmass

def ghost_plotting(choice, plot_size = (18, 6), plot_ratio = [10, 1], hist_bins = [50, 50, 10], labels = True, raw = False):
    if 'atlas_all' in choice:
        choice = ['altas_muvmass', 'altas_muvz', 'altas_muvmu', 'altas_residualsvz', 'altas_residualsvmass']
    if 'burns_all' in choice:
        choice = ['burns_muvmass', 'burns_muvz', 'burns_muvmu', 'burns_residualsvz', 'burns_residualsvmass']

    if ('burns_muvmass' in choice) or ('altas_muvmass' in choice):
        try:
            fig, axs = plt.subplots(1, 2, figsize=plot_size, gridspec_kw={'width_ratios': plot_ratio})

            if 'burns_muvmass' in choice:
                fig.suptitle('Host Mass of CSP 91bg-like SNe Ia') # Figure Title
                objs = dict_unpacker(BURNS_SAVE_TXT)
            elif 'altas_muvmass' in choice:
                fig.suptitle('Host Mass of ATLAS 91bg-like SNe Ia') # Figure Title
                objs = dict_unpacker(ATLAS_SAVE_TXT)

            mass_hist = []
            for obj in objs:
                mu, mu_err, mass = float(objs[obj]['mu']), float(objs[obj]['mu_err']), float(objs[obj]['logstellarmass'])
                if mass > 0:
                    axs[0].errorbar(mu, mass, xerr=mu_err, fmt='o')
                    if labels:
                        axs[0].text(mu, mass, obj, size='x-small', va='top')
                    mass_hist.append(mass)
            axs[1].hist(mass_hist, bins=hist_bins[0], orientation="horizontal")

            axs[0].set(xlabel='Distance Modulus', ylabel='Host Mass') # Sub-plot Labels
            axs[1].get_yaxis().set_visible(False)
            plt.tight_layout()
            plt.show()
        except KeyError:
            print("[KeyError] Please run 'ghost_host_galaxy()' before attempting to plot.")

    if ('burns_muvz' in choice) or ('altas_muvz' in choice):
        plot_scale = 'linear'
        sigma = 3

        try:
            fig, axs = plt.subplots(1, 2, figsize=plot_size, gridspec_kw={'width_ratios': plot_ratio})

            if 'burns_muvz' in choice:
                fig.suptitle('Distance Modulus vs. Redshift of CSP 91bg-like SNe Ia)') # Figure Title
                objs = dict_unpacker(BURNS_SAVE_TXT)
            elif 'altas_muvz' in choice:
                fig.suptitle('Distance Modulus vs. Redshift of ATLAS 91bg-like SNe Ia') # Figure Title
                objs = dict_unpacker(ATLAS_SAVE_TXT)

            # Plot -- mu vs z
            mu_hist = []
            for obj in objs:
                mu, mu_err, z = float(objs[obj]['mu']), float(objs[obj]['mu_err']), float(objs[obj]['z'])
                if not raw and mu <= 0:
                    continue
                axs[0].errorbar(z, mu, yerr=mu_err*sigma, fmt='o')
                if labels:
                    axs[0].text(z, mu, obj, size='x-small', va='top')
                mu_hist.append(mu)
            axs[1].hist(mu_hist, bins=hist_bins[1], orientation="horizontal")

            # Plot cosmology -- Use the distmod method of the cosmology object
            z_cosmo = np.linspace(0.001, 0.08, 100)
            axs[0].plot(z_cosmo, COSMO_MODEL.distmod(z_cosmo), alpha=0.4, linewidth=5,
                        label='H0:'+str(COSMO_MODEL._H0)+'\nOm0: '+str(COSMO_MODEL._Om0)+'\nOde0: '+str(COSMO_MODEL._Ode0))

            # Formating
            axs[0].set(xlabel='Redshift', ylabel='Distance Modulus') # Sub-plot Labels
            axs[0].set_xscale(plot_scale); axs[0].set_yscale(plot_scale)
            axs[0].legend()
            axs[1].get_yaxis().set_visible(False)
            fig.tight_layout()

            plt.show()
        except KeyError:
            print("[KeyError] Please run 'snpy_fit()' before attempting to plot.")

    if ('burns_muvmu' in choice) or ('altas_muvmu' in choice):
        plot_scale = 'linear'
        plt.figure(figsize=plot_size)

        try:
            if 'burns_muvmu' in choice:
                plt.title('SNooPy Distance Modulus vs. Cosmological Distance Modulus of CSP 91bg-like SNe Ia)') # Figure Title
                objs = dict_unpacker(BURNS_SAVE_TXT)
            elif 'altas_muvmu' in choice:
                plt.title('SNooPy Distance Modulus vs. Cosmological Distance Modulus of ATLAS 91bg-like SNe Ia') # Figure Title
                objs = dict_unpacker(ATLAS_SAVE_TXT)

            # Plot -- mu vs mu
            x_mu, y_mu = [], []
            for obj in objs:
                mu_snpy, z = float(objs[obj]['mu']), float(objs[obj]['z'])
                mu_cosmo = COSMO_MODEL.distmod(z).value
                plt.scatter(mu_snpy, mu_cosmo, marker='o')
                x_mu.append(mu_snpy)
                y_mu.append(mu_cosmo)
                if labels:
                    axs[0].text(mu_snpy, mu_cosmo, obj, size='x-small', va='top')

            # Plot trendline
            z = np.polyfit(x_mu, y_mu, 1)
            p = np.poly1d(z)
            plt.plot(x_mu, p(x_mu), alpha=0.4, linewidth=5, label='Trendline: '+"y=%.6fx+(%.6f)"%(z[0],z[1]))

            # Formating
            plt.xlabel('SNooPy Distance Modulus'); plt.ylabel('Cosmological Distance Modulus')
            plt.legend()

            plt.show()
        except KeyError:
            print("[KeyError] Please run 'snpy_fit()' before attempting to plot.")

    if ('burns_residualsvz' in choice) or ('altas_residualsvz' in choice):
        sigma = 0

        try:
            fig, axs = plt.subplots(1, 2, figsize=plot_size, gridspec_kw={'width_ratios': plot_ratio})
            if 'burns_residualsvz' in choice:
                fig.suptitle('Hubble Residuals vs. Redshift of CSP 91bg-like SNe Ia\n '+
                             'Dist. Sigma: '+str(sigma)) # Figure Title
                objs = dict_unpacker(BURNS_SAVE_TXT)
            elif 'altas_residualsvz' in choice:
                fig.suptitle('Hubble Residuals vs. Redshift of ATLAS 91bg-like SNe Ia\n'+
                             'Dist. Sigma: '+str(sigma)) # Figure Title
                objs = dict_unpacker(ATLAS_SAVE_TXT)

            # Compute mu_cosmo
            mu_cosmo, mu_snpy, mu_err, z, objnames = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
            for obj in objs:
                objnames = np.append(objnames, obj)
                z = np.append(z, float(objs[obj]['z']))
                mu_err = np.append(mu_err, float(objs[obj]['mu_err']))
                mu_snpy = np.append(mu_snpy, float(objs[obj]['mu']))
                mu_cosmo = np.append(mu_cosmo, COSMO_MODEL.distmod(float(objs[obj]['z'])).value)
            mu_res = (mu_snpy - mu_cosmo) - np.average(mu_snpy - mu_cosmo)

            # Plot -- mu vs z
            mu_hist = []
            for n in range(len(objnames)):
                if not raw and z[n] < 0.01:
                    continue
                axs[0].errorbar(z[n], mu_res[n], yerr=mu_err[n]*sigma, label=objnames[n], fmt="o")
                mu_hist.append(mu_res[n])
                if labels:
                    axs[0].text(z[n], mu_res[n], objnames[n], size='x-small', va='top')
            axs[1].hist(mu_hist, bins=hist_bins[2], orientation="horizontal")

            # Formatting
            axs[0].set(xlabel='Redshift', ylabel='Hubble Residuals') # Sub-plot Labels
            axs[1].get_yaxis().set_visible(False)
            plt.tight_layout()

            # Limits
            ylimiter = (np.max(np.abs(mu_hist)) + np.max(mu_err*sigma)) + 0.01
            axs[0].set_ylim(-ylimiter, ylimiter); axs[1].set_ylim(-ylimiter, ylimiter)

            plt.show()
        except KeyError:
            print("[KeyError] Please run 'snpy_fit()' before attempting to plot.")

    if ('burns_residualsvmass' in choice) or ('altas_residualsvmass' in choice):
        sigmas = [3, 3]
        objsRemove = ['2022skw', '2023cvq']

        try:
            fig, axs = plt.subplots(1, 2, figsize=plot_size, gridspec_kw={'width_ratios': plot_ratio})

            if 'burns_residualsvmass' in choice:
                fig.suptitle('Hubble Residuals vs. Host Mass of CSP 91bg-like SNe Ia\n '+
                             'Dist. Sigma: '+str(sigmas[0])+' | Mass Sigma: '+str(sigmas[1])) # Figure Title
                objs = dict_unpacker(BURNS_SAVE_TXT)
            elif 'altas_residualsvmass' in choice:
                fig.suptitle('Hubble Residuals vs. Host Mass of ATLAS 91bg-like SNe Ia\n '+
                             'Dist. Sigma: '+str(sigmas[0])+' | Mass Sigma: '+str(sigmas[1])) # Figure Title
                objs = dict_unpacker(ATLAS_SAVE_TXT)

            # Compute mu_cosmo
            mu_cosmo, mu_snpy, mu_err, mass, mass_err, objnames = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
            for obj in objs:
                if obj in objsRemove:
                    continue
                objnames = np.append(objnames, obj)
                mass = np.append(mass, float(objs[obj]['logstellarmass']))
                mass_err = np.append(mass_err, float(objs[obj]['logstellarmasserr']))
                mu_err = np.append(mu_err, float(objs[obj]['mu_err']))
                mu_snpy = np.append(mu_snpy, float(objs[obj]['mu']))
                mu_cosmo = np.append(mu_cosmo, COSMO_MODEL.distmod(float(objs[obj]['z'])).value)
            mu_res = (mu_snpy - mu_cosmo) - np.average(mu_snpy - mu_cosmo)

            # Plot -- mu vs z
            mu_hist = []
            for n in range(len(objnames)):
                if not raw and mass[n] <= 0:
                    continue
                axs[0].errorbar(mass[n], mu_res[n], xerr=mass_err[n]*sigmas[1], yerr=mu_err[n]*sigmas[0], fmt='o')
                mu_hist.append(mu_res[n])
                if labels:
                    axs[0].text(mass[n], mu_res[n], objnames[n], size='x-small', va='top')
            axs[1].hist(mu_hist, bins=hist_bins[3], orientation="horizontal")

            # Formatting
            axs[0].set(xlabel='Host Mass', ylabel='Hubble Residuals') # Sub-plot Labels
            axs[1].get_yaxis().set_visible(False)
            plt.tight_layout()

            # Limits
            ylimiter = (np.max(np.abs(mu_hist)) + np.max(mu_err*sigmas[0])) + 0.01
            axs[0].set_ylim(-ylimiter, ylimiter); axs[1].set_ylim(-ylimiter, ylimiter)
            # axs[0].set_xlim(9.75, 11.5); axs[1].set_xlim(9.75, 11.5)

            plt.show()
        except KeyError:
            print("[KeyError] Please run 'snpy_fit()' before attempting to plot.")

    if ('burns_res_zcorr' in choice) or ('altas_res_zcorr' in choice):
        sigma = 0
        objsRemove = ['2022skw', '2023cvq']

        try:
            fig, axs = plt.subplots(1, 2, figsize=plot_size, gridspec_kw={'width_ratios': plot_ratio})
            if 'burns_res_zcorr' in choice:
                fig.suptitle('Hubble Residuals vs. Redshift of CSP 91bg-like SNe Ia\n '+
                             'Dist. Sigma: '+str(sigma)) # Figure Title
                objs = dict_unpacker(BURNS_SAVE_TXT)
            elif 'altas_res_zcorr' in choice:
                fig.suptitle('Hubble Residuals vs. Redshift of ATLAS 91bg-like SNe Ia\n'+
                             'Dist. Sigma: '+str(sigma)) # Figure Title
                objs = dict_unpacker(ATLAS_SAVE_TXT)

            # Compute mu_cosmo
            mu_cosmo, mu_snpy, mu_err, z, objnames = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
            for obj in objs:
                if not raw and float(objs[obj]['z']) < 0.01:
                    continue
                elif obj in objsRemove:
                    continue
                objnames = np.append(objnames, obj)
                z = np.append(z, float(objs[obj]['z']))
                mu_err = np.append(mu_err, float(objs[obj]['mu_err']))
                mu_snpy = np.append(mu_snpy, float(objs[obj]['mu']))
                mu_cosmo = np.append(mu_cosmo, COSMO_MODEL.distmod(float(objs[obj]['z'])).value)
            mu_res = (mu_snpy - mu_cosmo) - np.average(mu_snpy - mu_cosmo)

            # Plot -- mu vs z
            mu_hist = []
            for n in range(len(objnames)):

                axs[0].errorbar(z[n], mu_res[n], yerr=mu_err[n]*sigma, label=objnames[n], fmt="o")
                mu_hist.append(mu_res[n])
                if labels:
                    axs[0].text(z[n], mu_res[n], objnames[n], size='x-small', va='top')
            axs[1].hist(mu_hist, bins=hist_bins[2], orientation="horizontal")

            # Formatting
            axs[0].set(xlabel='Redshift', ylabel='Hubble Residuals') # Sub-plot Labels
            axs[1].get_yaxis().set_visible(False)
            plt.tight_layout()

            # Limits
            ylimiter = (np.max(np.abs(mu_hist)) + np.max(mu_err*sigma)) + 0.1
            axs[0].set_ylim(-ylimiter, ylimiter); axs[1].set_ylim(-ylimiter, ylimiter)

            # Redshift correction
            # chi2 = malmquist_bias_corr(mu_res, z, mu_err, m, b)






            plt.close()
            # plt.show()
        except KeyError:
            print("[KeyError] Please run 'snpy_fit()' before attempting to plot.")


    return

""" GENERAL """

def recover_dir():
    directories = [ROOT_PATH,
                   DATA_ATLAS,
                   SNPY_ROOT, SNPY_BURNS, SNPY_BURNS_PLOTS,
                   SNPY_ATLAS, SNPY_ATLAS_PLOTS, SNPY_ATLAS_ASCII,
                   PLOTS_ROOT, TEST_ROOT, GHOST_DATA]
    files = [PROB_CHILD_TXT, TNS_KEY_TXT,
             BURNS_SAVE_TXT, ATLAS_SAVE_TXT]
    for dir in directories:
        if os.path.exists(dir) == False:
            os.mkdir(dir)
    for n_file in files:
        if os.path.exists(n_file) == False:
            with open(n_file, 'w') as f:
                pass

    astro_ghost.ghostHelperFunctions.getGHOST(real=False, verbose=True, installpath=GHOST_DATA, clobber=False)

    return


def TNS_details(ra, dec):
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


def snpy_fit(path, objname, save_loc, plot_save_loc, skip_problems=False, use_saved=False, snpy_plots=True,
             save_plots=True):
    problem_children = handle_problem_children(state='READ')  # Open problem children

    if skip_problems and (objname in problem_children):
        return None
    else:
        try:
            if use_saved and os.path.exists(save_loc + objname + '_EBV_model2.snpy'):
                n_s = snpy.get_sn(save_loc + objname + '_EBV_model2.snpy')
            else:
                n_s = snpy.get_sn(path)
                n_s.choose_model('EBV_model2', stype='st')
                n_s.set_restbands()  # Auto pick appropriate rest-bands

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
                print('\t', filter_wheel)

                n_s.fit(bands=filter_wheel, dokcorr=True, k_stretch=False, reset_kcorrs=True,
                        **{'mangle': 1, 'calibration': 0})
                n_s.save(save_loc + objname + '_EBV_model2.snpy')

            if snpy_plots:
                n_s.plot(outfile=plot_save_loc + objname + '_snpyplots.png')
                plt.show()
        except Exception as error:
            problem_children = np.append(problem_children, objname)
            handle_problem_children(state='WRITE', problem_c=problem_children)  # Commit problem children
            return error

    plt.close()
    return {'ra': n_s.ra, 'dec': n_s.decl, 'z': n_s.z, 'MJDs': min(mjds), 'MJDe': max(mjde),
            'mu': n_s.parameters['DM'], 'st': n_s.parameters['st'], 'Tmax': n_s.parameters['Tmax'],
            'EBVhost': n_s.parameters['EBVhost'],
            'mu_err': n_s.errors['DM'], 'st_err': n_s.errors['st'], 'Tmax_err': n_s.errors['Tmax'],
            'EBVhost_err': n_s.errors['EBVhost']}
def snpy_fit_indv(objname):
    path = '/content/HiloCATsSN1991bg/snpy/atlas/' + objname + '_EBV_model2.snpy'

    n_s = snpy.get_sn(path)
    n_s.choose_model('EBV_model2', stype='st')
    n_s.set_restbands()  # Auto pick appropriate rest-bands

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

    n_s.fit(bands=filter_wheel, dokcorr=True, k_stretch=False, reset_kcorrs=True, **{'mangle': 1, 'calibration': 0})
    plt.show()

    print('Results:',
          'mu =', n_s.parameters['DM'], '+/-', n_s.errors['DM'], '\n',
          '\t st =', n_s.parameters['st'], '+/-', n_s.errors['st'], '\n',
          '\t Tmax =', n_s.parameters['Tmax'], '+/-', n_s.errors['Tmax'], '\n',
          '\t EBVhost =', n_s.parameters['EBVhost'], '+/-', n_s.errors['EBVhost'], '\n',
          '\t MJD min:', min(mjds), '| MJD max:', max(mjde))

    return


def write_ASCII(objs, save_loc, quiet=True):
    print('[+++] Saving data to ASCII files for SNooPy...')
    for obj in objs:
        with open(save_loc + obj + '_snpy.txt', 'w') as f:
            # Line 1 -- Objname, Helio-Z, RA, Dec (Ex. SN1981D 0.005871 50.65992 -37.23272)
            f.write(
                str(obj) + ' ' + str(objs[obj]['z']) + ' ' + str(objs[obj]['ra']) + ' ' + str(objs[obj]['dec']) + '\n')

            # 'o'/'ATri'-filter photometry block -- Date (JD/MJD), mag, err (674.8593 12.94 0.11)
            f.write('filter ATri\n')
            for i in range(len(objs[obj]['time_o'])):
                f.write(str(objs[obj]['time_o'][i]) + '\t' + str(objs[obj]['mag_o'][i]) + '\t' + str(
                    objs[obj]['dmag_o'][i]) + '\n')

            # # 'c'/'ATgr'-filter photometry block
            f.write('filter ATgr\n')
            for i in range(len(objs[obj]['time_c'])):
                f.write(str(objs[obj]['time_c'][i]) + '\t' + str(objs[obj]['mag_c'][i]) + '\t' + str(
                    objs[obj]['dmag_c'][i]) + '\n')
    return


def read_DR3(loc='/content/HiloCATsSN1991bg/DR3_fits.dat'):
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
    catagories = list(data_dict[list(data_dict.keys())[0]].keys())
    with open(save_loc, 'w') as f:
        f.write('objname')
        for category in catagories:
            f.write(delimiter + category)
        f.write('\n')
        for objname in data_dict:
            f.write(objname)
            for category in catagories:
                f.write(delimiter + str(data_dict[objname][category]))
            f.write('\n')
    return


def handle_problem_children(state, problem_c=None):
    if state == 'READ':
        # Read problem children
        problem_c = np.genfromtxt(PROB_CHILD_TXT, dtype=str)
        return problem_c
    elif state == 'WRITE':
        # Write problem children
        problem_c = np.unique(problem_c)
        with open(PROB_CHILD_TXT, 'w') as f:
            for c in problem_c:
                f.write(c + '\n')
        return None
    else:
        raise Exception("Invalid state: '" + state + "' [READ/WRITE]")


def save_to_zip(zip_loc, save_loc):
    print('Saving zipped files to...', save_loc)
    files = glob.glob(zip_loc + '*')
    with ZipFile(save_loc, 'w') as zip:
        for n_file in files:
            zip.write(n_file)
    return


def malmquist_bias_corr(mu, z, mu_err, m, b):
    summation = 0
    for n in range(len(mu)):
        summation = summation + (mu[n] - ((m * z[n]) + b)) / (mu_err[n] ** 2)
    return summation


""" MAIN """
if __name__ == '__main__':
    start = systime.time() # Runtime tracker

    recover_dir() # Recovering vital directories

    # burns_cspvdr3(skip_problems=False, use_saved=False, snpy_plots=False)
    # ghost_host_galaxy(BURNS_SAVE_TXT, save_loc=TEST_ROOT, keep_data=False, update_saved=True)

    # atlas_collection(quiet=False, check_data=True)
    # atlas_objs = atlas_processing(err_max=1000, n_iter=0, sleep_t=5, use_TNS=True)
    # write_ASCII(atlas_objs, SNPY_ATLAS_ASCII, quiet=True)
    # atlas_snpy_fitting(n_iter=0, skip_problems=False, use_saved=False, snpy_plots=True, save_plots=True)
    # ghost_host_galaxy(ATLAS_SAVE_TXT, save_loc=TEST_ROOT, keep_data=True, update_saved=True)

    # ztf_collection(submit=True)
    # ztf_alt_collection()

    # burns_plotting([]) # Options: ['reg_hist', 'res_hist', 'zvmu']

    # atlas_plotting([])

    ghost_plot_args = {'plot_size': (18, 6), 'plot_ratio': [10, 1], 'hist_bins': [40, 40, 35, 35], 'labels': False, 'raw': False}
    ghost_plotting(['altas_residualsvz'], **ghost_plot_args)

    # snpy_fit_indv('2023cvq')

    print('|---------------------------|\n Run-time: ', round(systime.time()-start, 4), 'seconds\n|---------------------------|')