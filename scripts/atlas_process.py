import os
import glob
import numpy as np
import requests
import urllib.request
import scripts.general as gen
import time as systime
import matplotlib.pyplot as plt

DATA_ATLAS = '../data/ATLAS/'
SNPY_ATLAS = '../snpy/atlas/'
SNPY_ATLAS_PLOTS = SNPY_ATLAS+'plots/'
SNPY_ATLAS_ASCII = SNPY_ATLAS+'ascii/'
ATLAS_SAVE_TXT = SNPY_ATLAS+'atlas_saved.txt'

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
def atlas_processing(err_max=100, n_iter=10, sleep_t=5, use_TNS=False, loc_TNS='../working_data/TNS_key.txt'):
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
        tns_key = gen.dict_unpacker(loc_TNS)
        if use_TNS and (ATLAS_name in tns_key):
            print('\t\tAttempting to use TNS key for objname...')
            objname = tns_key[ATLAS_name]['objname']
            z = tns_key[ATLAS_name]['z']
        else:
            print('\t\tRetrieving TNS data for...', ATLAS_name, '[sleeping for', sleep_t, 'seconds...]')
            systime.sleep(sleep_t)
            details = gen.TNS_details(ra, dec)
            objname = details['objname']
            z = details['redshift']
            with open(loc_TNS, 'a') as f:
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
        temp_dict = gen.snpy_fit(objpaths[n], objname, save_loc=SNPY_ATLAS, plot_save_loc=SNPY_ATLAS_PLOTS, **fit_args)

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
        gen.dict_packer(objs, ATLAS_SAVE_TXT, delimiter=', ') # Save data from fitting
    if snpy_plots and save_plots:
        gen.save_to_zip(SNPY_ATLAS_PLOTS, SNPY_ATLAS_PLOTS+'atlas_snpy_plots.zip')

    return
def atlas_plotting(choice):
    if 'reg_hist' in choice:
        print('[+++] Ploting Histogram of ATLAS SNooPy Fitting Parameters...')
        objs = gen.dict_unpacker(ATLAS_SAVE_TXT)
        fig, ax = plt.subplots(1, 3, figsize=(15, 5))
        titles = ['st', 'ATLAS SNooPy Parameters' + '\nTmax', 'EBVhost']
        all_st, all_tmax, all_ebvhost = [], [], []
        for obj in objs:
            all_st.append(float(objs[obj]['st']))
            all_tmax.append(float(objs[obj]['Tmax']))
            all_ebvhost.append(float(objs[obj]['EBVhost']))
        ax[0].hist(all_st, bins=50, color='red')
        ax[1].hist(all_tmax, bins=50, color='purple')
        ax[2].hist(all_ebvhost, bins=50, color='blue')

        ax[0].set_title('st'); ax[1].set_title('ATLAS SNooPy Parameters' + '\nTmax'); ax[2].set_title('EBVhost')
        ax[1].set(xlabel='MJD+53000')

        plt.savefig('../save/plots/atlas_param_hist.png')

        plt.show()