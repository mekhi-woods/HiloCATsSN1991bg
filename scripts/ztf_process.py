import os
import numpy as np
import scripts.general as gen
import requests
import json
import matplotlib.pyplot as plt
import glob
import time as systime

DATA_ZTF = '../data/ZTF/'
SNPY_ZTF = '../snpy/ztf/'
SNPY_ZTF_PLOTS = SNPY_ZTF+'plots/'
SNPY_ZTF_ASCII = SNPY_ZTF+'ascii/'
ZTF_SAVE_TXT = SNPY_ZTF+'ztf_saved.txt'

def ztf_wget(submit=False):
    atlas_txt = np.genfromtxt('../txts/91bg_list.txt', dtype=str, skip_header=1)

    print('Requesting the following from ZTF...')
    for i in range(len(atlas_txt[:, 0])):
        print('---------------------------------------')
        ra = atlas_txt[i, 1]
        dec = atlas_txt[i, 2]

        # ra = atlas_txt[tar]['ra']
        # dec = atlas_txt[tar]['dec']
        # jds = float(atlas_txt[tar]['MJDs']) + 2400000 # MJD to JD
        # jde = float(atlas_txt[tar]['MJDe']) + 2400000

        if submit:
            # cmd = f"wget --http-user=ztffps --http-passwd=dontgocrazy! -O log.txt \"https://ztfweb.ipac.caltech.edu/cgi-bin/requestForcedPhotometry.cgi?ra={ra}&dec={dec}&jdstart={jds}&jdend={jde}&email={'mekhidw@hawaii.edu'}&userpass={'wxdk286'}\""
            cmd = f"wget --http-user=ztffps --http-passwd=dontgocrazy! -O log.txt \"https://ztfweb.ipac.caltech.edu/cgi-bin/requestForcedPhotometry.cgi?ra={ra}&dec={dec}&email={'mekhidw@hawaii.edu'}&userpass={'wxdk286'}\""
            # print(cmd)
            results = os.system(cmd)

            print('RA: '+str(ra)+', DEC: '+str(dec)+', '+str(results))

    return
def ztf_alt_processing(paths=glob.glob('../data/ZTF/old/*.txt'), min_pts=10, mag_err_max=0.75, flux_err_max=80):
    print('[+++] Retrieving data from...', DATA_ZTF)
    header_num = 57
    ZTFobjs = {}

    for path in paths:
        # Tracking/Cosmetics
        ztfname = path.split('/')[-1][11:-7].split('_')[-1][3:]
        tracker = '['+str(paths.index(path)+1)+'/'+str(len(paths))+']' # Purely cosmetic
        print(tracker, '\n', '\t\tPulling', ztfname, 'data...')

        data = np.genfromtxt(path, delimiter=' ', skip_header=header_num, dtype=str)

        # Grab header
        hdr = []
        with open(path, 'r') as f:
            for i in range(header_num):
                hdr.append(f.readline())

        # Get objname from TNS
        name_found = False
        tns_key = np.genfromtxt('../working_data/TNS_key.txt', dtype=str, delimiter=', ', skip_header=1)
        if len(tns_key) > 0:
            objnames, ra, dec, z = tns_key[:, 1], tns_key[:, 2].astype(float), tns_key[:, 3].astype(float), tns_key[:,4]
            obj_ra, obj_dec = float(hdr[3].split(' ')[5]), float(hdr[4].split(' ')[5])
            for n in range(len(objnames)):
                if abs(obj_ra - ra[n]) < 0.01:
                    obj_name, obj_z = objnames[n], z[n]
                    name_found = True
                    break
        if name_found == False:
            details = gen.TNS_details(hdr[3].split(' ')[5], hdr[4].split(' ')[5])
            obj_name, obj_ra, obj_dec, obj_z = details['objname'], details['radeg'], details['decdeg'], details[
                'redshift']
            with open('../working_data/TNS_key.txt', 'a') as f:
                f.write(ztfname + ', ' + obj_name + ', ' + str(obj_ra) + ', ' + str(obj_dec) + ', ' + str(obj_z))

        obj_filters, obj_time, obj_flux, obj_flux_unc, obj_mag, obj_mag_unc = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
        for i in range(len(data[:, 0])):
            # Remove intial nulls
            valid = True
            for j in [4, 22, 24, 25]:
                if data[i, j] == 'null':
                    valid = False
            if not valid:
                continue

            # Flux to magnitude
            n_flux, n_flux_unc = float(data[i, 24]), float(data[i, 25])
            n_mag = -2.5 * np.log10(n_flux) + 23.9
            n_mag_unc = (2.5 / (n_flux * np.log(10))) * n_flux_unc
            if np.isnan(n_mag):
                continue

            # Clean data
            if n_flux <= 0: # Negatives
                continue
            if n_flux_unc > flux_err_max: # Extreme Flux
                continue
            if n_mag_unc > mag_err_max: # Extreme Mags
                continue

            # Commit variables
            # print(i, n_flux, n_flux_unc, n_mag, n_mag_unc)
            obj_filters = np.append(obj_filters, data[i, 4])
            obj_time = np.append(obj_time, float(data[i, 22]))
            obj_flux = np.append(obj_flux, n_flux)
            obj_flux_unc = np.append(obj_flux_unc, n_flux_unc)
            obj_mag = np.append(obj_mag, n_mag)
            obj_mag_unc = np.append(obj_mag_unc, n_mag_unc)

        # Create Object
        # flux(forcediffimflux), flux_unc(forcediffimfluxunc)
        obj = {'ra': obj_ra, 'dec': obj_dec, 'z': obj_z, 'filters': obj_filters, 'time': obj_time,
               'flux': obj_flux, 'flux_unc': obj_flux_unc, 'mag': obj_mag, 'mag_unc': obj_mag_unc}
        print('\t\tRetrived: ' + obj_name + ' | ra: '+str(obj_ra)+' \ dec: '+str(obj_dec))

        # Commit Object
        ZTFobjs.update({obj_name: obj})

    print('[!!!]\t\tRetrieved & processed', len(ZTFobjs), 'SNe from ZTF!')

    return ZTFobjs
def ztf_plotting(ZTFobjs, choice = 'flux', plot_filters=['g', 'r', 'i'], zoom = 330, sigma = 1, pause_time=3):
    i = 0
    for obj in ZTFobjs:
        i += 1
        tracker = '['+str(i)+' / '+str(len(ZTFobjs))+']'
        print(tracker, 'Ploting '+obj+'...')
        plt.figure(figsize=(12, 6))

        # Filter selection
        g_filters = np.where(ZTFobjs[obj]['filters'] == 'ZTF_g')[0]
        r_filters = np.where(ZTFobjs[obj]['filters'] == 'ZTF_r')[0]
        i_filters = np.where(ZTFobjs[obj]['filters'] == 'ZTF_i')[0]

        if choice == 'flux':
            if 'g' in plot_filters:
                plt.errorbar(ZTFobjs[obj]['time'][g_filters], ZTFobjs[obj]['flux'][g_filters], yerr=ZTFobjs[obj]['flux_unc'][g_filters] * sigma, fmt='o',
                             color='g', label='g-fitler')
            if 'r' in plot_filters:
                plt.errorbar(ZTFobjs[obj]['time'][r_filters], ZTFobjs[obj]['flux'][r_filters], yerr=ZTFobjs[obj]['flux_unc'][r_filters] * sigma, fmt='o',
                             color='r', label='r-fitler')
            if 'i' in plot_filters:
                plt.errorbar(ZTFobjs[obj]['time'][i_filters], ZTFobjs[obj]['flux'][i_filters], yerr=ZTFobjs[obj]['flux_unc'][i_filters] * sigma, fmt='o',
                             color='violet', label='i-fitler')

            # Zooming in on Tmax
            if zoom > 0:
                f_max_index = np.where(ZTFobjs[obj]['flux'] == np.max(ZTFobjs[obj]['flux']))[0][0]
                t_max = ZTFobjs[obj]['time'][f_max_index]
                plt.xlim(t_max - zoom, t_max + zoom)
        if choice == 'mag':
            if 'g' in plot_filters:
                plt.errorbar(ZTFobjs[obj]['time'][g_filters], ZTFobjs[obj]['mag'][g_filters],
                             yerr=ZTFobjs[obj]['mag_unc'][g_filters] * sigma, fmt='o',
                             color='g', label='g-fitler')
            if 'r' in plot_filters:
                plt.errorbar(ZTFobjs[obj]['time'][r_filters], ZTFobjs[obj]['mag'][r_filters],
                             yerr=ZTFobjs[obj]['mag_unc'][r_filters] * sigma, fmt='o',
                             color='r', label='r-fitler')
            if 'i' in plot_filters:
                plt.errorbar(ZTFobjs[obj]['time'][i_filters], ZTFobjs[obj]['mag'][i_filters],
                             yerr=ZTFobjs[obj]['mag_unc'][i_filters] * sigma, fmt='o',
                             color='violet', label='i-fitler')

            # Zooming in on Tmax
            if zoom > 0:
                f_max_index = np.where(ZTFobjs[obj]['mag'] == np.min(ZTFobjs[obj]['mag']))[0][0]
                t_max = ZTFobjs[obj]['time'][f_max_index]
                plt.xlim(t_max - zoom, t_max + zoom)
            plt.gca().invert_yaxis()

        # Formating
        plt.xlabel('MJD+58000');
        plt.ylabel('Forced Difference Image Flux')
        plt.title(obj +
                  '\n(' + str(round(ZTFobjs[obj]['ra'], 5)) + ', ' + str(round(ZTFobjs[obj]['dec'], 5)) + ')' +
                  '\n Sigma: ' + str(sigma) + ', Zoom: ' + str(zoom) + ' days')
        plt.legend()

        plt.show()
        systime.sleep(pause_time)

    return
def ztf_write_ASCII(ZTFobjs, save_loc, quiet=True):
    print('[+++] Saving data to ASCII files for SNooPy...')
    for obj in ZTFobjs:
        with open(save_loc+obj+'_snpy.txt', 'w') as f:
            # Line 1 -- Objname, Helio-Z, RA, Dec (Ex. SN1981D 0.005871 50.65992 -37.23272)
            f.write(str(obj)+' '+str(ZTFobjs[obj]['z'])+' '+str(ZTFobjs[obj]['ra'])+' '+str(ZTFobjs[obj]['dec'])+'\n')

            # 'g'-filter photometry block -- Date (JD/MJD), mag, err (674.8593 12.94 0.11)
            g_filters = np.where(ZTFobjs[obj]['filters'] == 'ZTF_g')[0]
            f.write('filter g\n')
            for i in range(len(ZTFobjs[obj]['time'][g_filters])):
                f.write(str(ZTFobjs[obj]['time'][g_filters][i]+58000)+'\t'+str(ZTFobjs[obj]['mag'][g_filters][i])+'\t'+str(ZTFobjs[obj]['mag_unc'][g_filters][i])+'\n')

            # 'r'-filter photometry block
            r_filters = np.where(ZTFobjs[obj]['filters'] == 'ZTF_r')[0]
            f.write('filter r\n')
            for i in range(len(ZTFobjs[obj]['time'][r_filters])):
                f.write(str(ZTFobjs[obj]['time'][r_filters][i]+58000)+'\t'+str(ZTFobjs[obj]['mag'][r_filters][i])+'\t'+str(ZTFobjs[obj]['mag_unc'][r_filters][i])+'\n')

            # 'i'-filter photometry block
            i_filters = np.where(ZTFobjs[obj]['filters'] == 'ZTF_i')[0]
            f.write('filter i\n')
            for i in range(len(ZTFobjs[obj]['time'][i_filters])):
                f.write(str(ZTFobjs[obj]['time'][i_filters][i]+58000) + '\t' + str(ZTFobjs[obj]['mag'][i_filters][i]) + '\t' + str(ZTFobjs[obj]['mag_unc'][i_filters][i]) + '\n')

    return
def ztf_snpy_fitting(n_iter=0, skip_problems=True, use_saved=True, snpy_plots=True, save_plots=True):
    print('[+++] Fitting ZTF data with SNooPy...')
    fit_args = {'skip_problems': skip_problems, 'use_saved': use_saved, 'snpy_plots': snpy_plots, 'save_plots': save_plots}
    print('Fitting arguments: ', fit_args)
    objpaths = glob.glob(SNPY_ZTF_ASCII+'*')
    if n_iter != 0 and n_iter <= len(objpaths):
        objpaths = objpaths[:n_iter]

    objs = {}
    err_i = 0
    for n in range(len(objpaths)):

        # # rm after testing
        # if n+1 != 5:
        #     continue

        tracker = '['+str(n+1)+'/'+str(len(objpaths))+']' # Purely cosmetic
        objname = objpaths[n][len(SNPY_ZTF_ASCII):-9]
        print(tracker, objname)
        temp_dict = gen.newer_snpy_fit(objpaths[n], objname, save_loc=SNPY_ZTF, plot_save_loc=SNPY_ZTF_PLOTS, **fit_args)

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

    print('Finshed! Successfully fit', len(objpaths)-err_i, '/', len(objpaths), 'SNe from ZTF! ['+str(round(((len(objpaths)-err_i) / len(objpaths))*100, 2))+'%]')

    # Save Data
    if len(objs) > 0:
        gen.dict_packer(objs, ZTF_SAVE_TXT, delimiter=', ') # Save data from fitting
    if snpy_plots and save_plots:
        gen.save_to_zip(SNPY_ZTF_PLOTS, SNPY_ZTF_PLOTS+'ztf_snpy_plots.zip')

    return









