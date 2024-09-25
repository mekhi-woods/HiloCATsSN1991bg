import warnings
warnings.simplefilter("ignore", UserWarning)
import os
import sys
import snpy
import glob
import shutil
import sncosmo
import numpy as np
import time as systime
import matplotlib.pyplot as plt
from astro_ghost.ghostHelperFunctions import getTransientHosts
from astropy.coordinates import SkyCoord, Galactic
from astropy.table import Table
from astropy import units as u
from astroquery.sdss import SDSS
from scripts import general as gen

CONSTANTS = gen.get_constants()

class sn91bg():
    def __init__(self, objname=None, originalname=None, coords=(0.00, 0.00), z=0.00, origin=None, discovery_data=None):
        self.objname = objname
        self.originalname = originalname
        self.coords = coords
        self.z = z
        self.z_cmb = np.nan
        self.origin = origin
        self.discovery_date = discovery_data

        self.period = None
        self.params = {}

        self.zp = np.array([])
        self.filters = np.array([])
        self.time = np.array([])
        self.flux = np.array([])
        self.dflux = np.array([])
        self.mag = np.array([])
        self.dmag = np.array([])

        return
    def __str__(self):
        return self.origin + ' | ' + self.objname + ' | ' + self.originalname + ' | ' + str(self.coords) + ' | z = ' + str(self.z)
    def print_info(self):
        prnt_str = ('---------------------------------------------------------------------------------------------\n' +
                    self.origin + ' | ' + self.objname + ' | ' + str(self.coords) +
                    ' | z = ' + str(self.z) + ' | z_cmb = ' + str(self.z_cmb) + '\n' +
                    '---------------------------------------------------------------------------------------------\n' +
                    '\tFilters = ' + str(np.unique(self.filters)) + '\n' +
                    '\tFlux Range = (' + str(np.min(self.flux)) + ' +/- ' + str(self.dflux[np.where(self.flux == np.min(self.flux))[0][0]]) +
                    ') -- (' + str(np.max(self.flux)) + ' +/- ' + str(self.dflux[np.where(self.flux == np.max(self.flux))[0][0]]) + ')\n' +
                    '\tMagnitude Range = (' + str(np.min(self.mag)) + ' +/- ' + str(self.dmag[np.where(self.mag == np.min(self.mag))[0][0]]) +
                    ') -- (' + str(np.max(self.mag)) + ' +/- ' + str(self.dmag[np.where(self.mag == np.max(self.mag))[0][0]]) + ')\n' +
                    '\tMJDs-MJDe = ' + str(self.period[0]) + ' -- ' + str(self.period[1]) + '\n' +
                    '---------------------------------------------------------------------------------------------\n')
        if len(self.params) > 0:
            for p in self.params:
                prnt_str += ('\t'+p+' = ' + str(self.params[p]['value']) + ' +/- ' + str(self.params[p]['err']) + '\n')
            prnt_str += '---------------------------------------------------------------------------------------------'
        print(prnt_str)
    def plot(self, y_type='mag', save_loc='', zoom=0, subplots=False, date_lines=True):
        print('[+++] Plotting LC of '+self.objname+'...')
        filter_dict = {'u': 'teal', 'g': 'green', 'r': 'red', 'i': 'indigo', 'B': 'blue',
                       'V0': 'violet', 'V1': 'purple', 'V': 'red', 'Y': 'goldenrod', 'Hdw': 'tomato', 'H': 'salmon',
                       'J': 'aquamarine', 'Jrc2': 'cadetblue', 'Jdw': 'turquoise', 'Ydw': 'olive',
                       'c': 'cyan', 'o': 'orange', 'ZTF_g': 'green', 'ZTF_r': 'red', 'ZTF_i': 'indigo'}

        # Plot
        unique_filters, num_plts = np.unique(self.filters), len(np.unique(self.filters))
        if subplots == False:
            fig = plt.figure(figsize=(10, 4), constrained_layout=True)
            size = None
        elif self.origin == 'CSP':
            fig = plt.figure(figsize=(25, 8), constrained_layout=True)
            size = (2, 6)
        elif self.origin == 'ATLAS':
            fig = plt.figure(figsize=(12, 4), constrained_layout=True)
            size = (1, 2)
        elif self.origin == 'ZTF':
            fig = plt.figure(figsize=(25, 5), constrained_layout=True)
            size = (1, 3)
        elif self.origin == 'ATLAS-ZTF':
            fig = plt.figure(figsize=(25, 8), constrained_layout=True)
            size = (2, 5)
        else:
            raise ValueError('[!!!] Origin not valid!'
                             )

        # Plot for each filter
        for i in range(num_plts):
            if size != None:
                plt.subplot(size[0], size[1], i+1)
            indexs = np.where(self.filters == unique_filters[i])[0]
            if y_type == 'mag':
                plt.errorbar(self.time[indexs], self.mag[indexs], yerr=self.dmag[indexs], fmt='o', ms=4, elinewidth=0.3,
                             color=filter_dict[self.filters[indexs][0]], label=self.filters[indexs][0])
            elif y_type == 'flux':
                plt.errorbar(self.time[indexs], self.flux[indexs], yerr=self.dflux[indexs], fmt='o', ms=4, elinewidth=0.3,
                             color=filter_dict[self.filters[indexs][0]], label=self.filters[indexs][0])

            # Format
            if size != None and i > 0:
                plt.gca().get_yaxis().set_visible(False)  # Turn off y-axis labels

        if date_lines:
            # Tmax line
            if len(self.params) > 0:
                if 'Tmax' in list(self.params.keys()):
                    plt.axvline(x=self.params['Tmax']['value'],
                                color='maroon', ls='-.', label='Tmax', linewidth=3, alpha=0.5)
                elif 't0' in list(self.params.keys()):
                    plt.axvline(x=self.params['t0']['value'],
                                color='maroon', ls='-.', label='Tmax', linewidth=3, alpha=0.5)

            # Plot discovery date
            plt.axvline(x=float(self.discovery_date), color='peru', ls='--', label='Discovery Date', linewidth=3, alpha=0.5)

        if zoom > 0:
            if 'Tmax' in list(self.params.keys()):
                plt.xlim(self.params['Tmax']['value']-zoom, self.params['Tmax']['value']+zoom)
            elif 't0' in list(self.params.keys()):
                plt.xlim(self.params['t0']['value']-zoom, self.params['t0']['value']+zoom)

            plt.xlabel('MJD'); plt.legend()
        if y_type == 'mag':
            plt.gca().invert_yaxis()
        elif y_type == 'flux':
            plt.ylim(0)
        plt.suptitle('Lightcurve -- '+self.objname+' | '+self.originalname+' -- '+y_type)
        plt.show()
        return
    # -----------------------------------------------------------------------------------------------------------------
    def save_class(self, save_loc):
        print('[+++] '+self.objname+' -- Saving class to '+save_loc+'classes/'+self.objname+'_class.txt')
        with open(save_loc+'classes/'+self.objname+'_class.txt', 'w') as f:
            f.write(self.origin + ' ' + self.objname + ' ' + self.originalname +
                    ' ' + str(self.coords[0]) + ' ' + str(self.coords[1]) +
                    ' ' + str(self.z) + ' ' + str(self.z_cmb) + ' ' + str(self.discovery_date) + '\n')
            f.write('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')

            for p in self.params:
                f.write(p+', ' + str(self.params[p]['value']) + ', ' + str(self.params[p]['err']) + '\n')

            f.write('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
            cats = [self.zp, self.filters, self.time, self.flux, self.dflux, self.mag, self.dmag]
            cat_names = ['zp', 'filters', 'time', 'flux', 'dflux', 'mag', 'dmag']
            for i in range(len(cat_names)):
                f.write(cat_names[i])
                for j in range(len(self.zp)):
                    f.write(', ' + str(cats[i][j]))
                f.write('\n')
        return
    def load_from_file(self, file):
        print('[+++] Loading class from '+file)
        with open(file, 'r') as f:
            hdr = f.readline().split(' ')
            self.origin, self.objname, self.originalname, self.coords, self.z, self.z_cmb, self.discovery_date = (
                hdr[0], hdr[1], hdr[2], (float(hdr[3]), float(hdr[4])), float(hdr[5]), float(hdr[6]), float(hdr[7][:-1]))

            skip = f.readline()

            line = f.readline()
            while '+++' not in line:
                line = line.split(', ')
                self.params.update({line[0]: {'value': float(line[1]), 'err': float(line[2])}})
                line = f.readline()

            self.zp = np.array(f.readline().split(', ')[1:]).astype(float)
            self.filters = np.array(f.readline().split(', ')[1:]).astype(str)
            self.filters[-1] = self.filters[-1][:-1] # Removes the /n from the end of the line
            self.time = np.array(f.readline().split(', ')[1:]).astype(float)
            self.flux = np.array(f.readline().split(', ')[1:]).astype(float)
            self.dflux = np.array(f.readline().split(', ')[1:]).astype(float)
            self.mag = np.array(f.readline().split(', ')[1:]).astype(float)
            self.dmag = np.array(f.readline().split(', ')[1:]).astype(float)

            self.period = (np.min(self.time), np.max(self.time))
        return self
    # -----------------------------------------------------------------------------------------------------------------
    def fit(self, algo='snpy', save_loc=None):
        if save_loc == None:
            save_loc = CONSTANTS[algo+'_'+self.origin.lower() + '_saved_loc']

        # Choose alogrithm
        if algo == 'snpy':
            self.write_snpy_ASCII(save_loc= save_loc+ 'ascii/')
            self.snpy_fit(save_loc=save_loc)
        elif algo == 'salt':
            self.salt_fit(save_loc=save_loc)
        if self.params['mu']['value'] <= 0.00:
            return None

        # Analysis
        self.get_host_mass(use_key=True)

        # Save
        self.save_class(save_loc)
        return
    def make_class(self, data_set, path, dmag_max=0, dflux_max=0):
        if data_set == 'CSP':
            print('[+++] Creating class using CSP data...')
            print(path)
            # Header elements
            with open(path, 'r') as f:
                objname, z, ra, dec = f.readline().split(' ')
                objname, z, ra, dec = objname[2:], float(z), float(ra), float(dec[:-1])
            originalname = path.split('/')[-1].split('_')[0]

            # Query TNS for transient details
            objname, z_void, discdate = gen.TNS_details(ra, dec) # Voiding the redshift from TNS, its None for some reason

            # Make class
            tempSN = sn91bg(objname=objname,
                            originalname=originalname,
                            coords = (ra, dec),
                            z=z,
                            origin = 'CSP',
                            discovery_data=discdate)

            # Set arrays
            filter = ''
            with open(path, 'r') as f:
                hdr = f.readline()
                for line in f.readlines():
                    data_line = line[:-1].split(' ')
                    if len(data_line) == 2:
                        filter = data_line[1]
                    elif len(data_line) >= 3:
                        if len(data_line) == 4:
                            data_line = data_line[1:]
                        n_time, n_mag, n_dmag = float(data_line[0]), float(data_line[1]), float(data_line[2])
                        zp = float(CONSTANTS['csp_zpts_' + filter])
                        n_time = n_time + 53000  # JD to MJD

                        n_flux = 10 ** ((n_mag - zp) / -2.5)
                        n_dflux = np.abs(n_flux) * np.log(10) * ((1 / 2.5) * n_dmag)

                        tempSN.zp = np.append(tempSN.zp, zp)
                        tempSN.filters = np.append(tempSN.filters, filter)
                        tempSN.time = np.append(tempSN.time, n_time)
                        tempSN.flux = np.append(tempSN.flux, n_flux)
                        tempSN.dflux = np.append(tempSN.dflux, n_dflux)
                        tempSN.mag = np.append(tempSN.mag, n_mag)
                        tempSN.dmag = np.append(tempSN.dmag, n_dmag)
        elif data_set == 'ATLAS':
            print('[+++] Creating class using ATLAS data...')
            print(path)

            # Load data
            data = np.genfromtxt(path, delimiter=',', dtype=str, skip_header=1)

            # Check if  file empty
            if len(data) == 0:
                print('[!!!] File [' + path + '] empty!')
                return None

            # Query TNS for transient details
            ra, dec = np.average(data[:, 1].astype(float)), np.average(data[:, 2].astype(float))
            objname, z, discdate = gen.TNS_details(ra, dec)

            # Make class
            tempSN = sn91bg(objname=objname,
                            originalname=path.split('/')[-1].split('.')[0],
                            coords = (ra, dec),
                            z=np.nan if z == 'None' else float(z),
                            origin = 'ATLAS',
                            discovery_data=discdate)

            # Set arrays
            tempSN.zp = data[:, 7].astype(float)
            tempSN.filters = data[:, 6]
            tempSN.time = data[:, 8]
            tempSN.flux = data[:, 16]
            tempSN.dflux = data[:, 17]
            tempSN.mag = data[:, 3]
            tempSN.dmag = data[:, 4]
        elif data_set == 'ZTF':
            print('[+++] Creating class using ZTF data...')
            print(path)

            # Load data
            data = np.genfromtxt(path, delimiter=None, dtype=str, skip_header=56)

            # Check if file is empty
            if len(data) == 0:
                print('[!!!] File [' + path + '] empty!')
                return None

            with (open(path, 'r') as f):
                ztf_spread = float(CONSTANTS['ztf_spread'])

                hdr = f.readlines()
                ra, dec = float(hdr[3].split(' ')[-2]), float(hdr[4].split(' ')[-2])
                objname, z, discdate = gen.TNS_details(ra, dec)
                originalname = path.split('/')[-1].split('.')[0].split('_')[1]
                z = np.nan if z == 'None' else float(z)

                # Get magnitudes m = -2.5log(F) + zp
                time, flux, dflux, zp, filters = data[:, 22], data[:, 24], data[:, 25], data[:, 20], data[:, 4]
                valid_ints = np.unique(np.hstack((np.where(flux != 'null')[0], np.where(dflux != 'null')[0])))
                time, zp, filters = time[valid_ints].astype(float), zp[valid_ints].astype(float), filters[valid_ints]
                flux, dflux = flux[valid_ints].astype(float), dflux[valid_ints].astype(float)
                time = time - 2400000.5  # JD to MJD
                mag = (-2.5 * np.log10(flux)) + zp
                dmag = np.abs(-1.08573620476 * (dflux / flux))

                # Adjusting around tmax
                if ztf_spread != 0 and len(time) != 0:
                    t_max_guess = float(discdate)
                    zoom_indexes = np.where(time < t_max_guess + ztf_spread)
                    zoom_indexes = np.where(time[zoom_indexes] > t_max_guess - ztf_spread)

                    time = time[zoom_indexes]
                    flux = flux[zoom_indexes]
                    dflux = dflux[zoom_indexes]
                    mag = mag[zoom_indexes]
                    dmag = dmag[zoom_indexes]
                    zp = zp[zoom_indexes]
                    filters = filters[zoom_indexes]

            # Make class
            tempSN = sn91bg(objname=objname,
                            originalname=originalname,
                            coords = (ra, dec),
                            z=z,
                            origin = 'ZTF',
                            discovery_data=discdate)

            # Set arrays
            tempSN.zp = zp
            tempSN.filters = filters
            tempSN.time = time
            tempSN.flux = flux
            tempSN.dflux = dflux
            tempSN.mag = mag
            tempSN.dmag = dmag
        else:
            raise ValueError("Data set '" +
                             data_set + "' not recognized")

        # Clean data
        tempSN.clean_data(dmag_max, dflux_max)

        # Final stage
        if tempSN.period == None:
            print('[!!!] No valid points found in file!')
            return None
        else:
            print('      Class created successfully!')
            return tempSN
    def clean_data(self, dmag_max=0, dflux_max=0):
        print('[+++] '+self.objname+' -- Cleaning data...')

        # Adjust Maximums
        if dmag_max == 'median':
            dmag_max = np.median(self.dmag)
        if dflux_max == 'median':
            dflux_max = np.median(self.dflux)
        if dmag_max == 'average':
            dmag_max = np.median(self.dmag)
        if dflux_max == 'average':
            dflux_max = np.average(self.dflux)

        new_zp, new_filters, new_time, new_mag, new_dmag, new_flux, new_dflux = (
            np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]))
        for n in range(len(self.zp)):
            n_zp, n_filters, n_time  = self.zp[n], self.filters[n], self.time[n]
            n_mag, n_dmag, n_flux, n_dflux = self.mag[n], self.dmag[n], self.flux[n], self.dflux[n]
            n_mag = str(n_mag).replace('>', '')
            if n_mag == 'None' or n_dmag == 'None' or n_flux == 'None' or n_dflux == 'None':
                continue
            if float(n_dflux) == 0 or float(n_dmag) == 0 or float(n_mag) <= 0 or float(n_flux) <= 0:
                continue
            # Cut errors
            if (dmag_max != 0) and (float(n_dmag) > dmag_max):
                continue
            if (dflux_max != 0) and (float(n_dflux) > dflux_max):
                continue

            new_zp = np.append(new_zp, float(n_zp))
            new_filters = np.append(new_filters, n_filters)
            new_time = np.append(new_time, float(n_time))
            new_mag = np.append(new_mag, float(n_mag))
            new_dmag = np.append(new_dmag, float(n_dmag))
            new_flux = np.append(new_flux, float(n_flux))
            new_dflux = np.append(new_dflux, float(n_dflux))

        self.zp = np.copy(new_zp)
        self.filters = np.copy(new_filters)
        self.time = np.copy(new_time)
        self.mag = np.copy(new_mag)
        self.dmag = np.copy(new_dmag)
        self.flux = np.copy(new_flux)
        self.dflux = np.copy(new_dflux)

        if len(self.time) > 0:
            self.period = (np.min(self.time), np.max(self.time))

        return
    def write_snpy_ASCII(self, save_loc='default/'):
        print('[+++] '+self.objname+' -- Saving data to ASCII files for SNooPy...')
        filter_dict = {'o': 'ATri', 'c': 'ATgr',
                       'ZTF_g': 'g', 'ZTF_r': 'r', 'ZTF_i': 'i',
                       'B': 'B', 'H': 'H', 'J': 'J', 'Jrc2': 'Jrc2', 'V': 'V', 'V0': 'V0', 'Y': 'Y', 'Ydw': 'Ydw', 'g': 'g', 'i': 'i', 'r': 'r', 'u': 'u'}
        with open(save_loc + self.objname + '_snpy.txt', 'w') as f:
            # Line 1 -- Objname, Helio-Z, RA, Dec (Ex. SN1981D 0.005871 50.65992 -37.23272)
            f.write(str(self.objname)+' '+str(self.z)+' '+str(self.coords[0])+' '+str(self.coords[1])+'\n')
            for f_w in np.unique(self.filters):
                f_indexs = np.where(self.filters == f_w)[0]
                f.write('filter ' + filter_dict[f_w] + '\n')
                for i in f_indexs:
                    # filter photometry block -- Date (JD/MJD), mag, err (i.e. 674.8593 12.94 0.11)
                    f.write(str(self.time[i]) + '\t' + str(self.mag[i]) + '\t' + str(self.dmag[i]) + '\n')
        print('      Saved file to '+save_loc + self.objname + '_snpy.txt')
        return
    def snpy_fit(self, save_loc, use_saved=True, show_plot=True, quiet=False):
        print('[+++] '+self.objname+' -- Fitting data with SNooPy...')
        load_path = save_loc + 'ascii/' + self.objname + '_snpy.txt'
        save_path = save_loc + 'models/' + self.objname + '_EBV_model2.snpy'
        param_names = ['mu', 'st', 'Tmax', 'EBVhost']
        snpy_param_names = ['DM', 'st', 'Tmax', 'EBVhost']

        # Check quiet
        if quiet:
            sys.stdout = open(os.devnull, 'w')

        # Check saved models
        if use_saved and os.path.isfile(save_path):
            print('[+++] Saved model found! Pulling from...', save_path)
            n_s = snpy.get_sn(save_path)
            for i in range(len(param_names)):
                self.params.update({param_names[i]: {'value': n_s.parameters[snpy_param_names[i]],
                                                     'err': n_s.errors[snpy_param_names[i]]}})
            return

        # Load Data
        try:
            n_s = snpy.get_sn(load_path)
        except:
            self.params.update({'mu': {'value': 0.00, 'err': 0.00}})
            print('[!!!] Failed to load ASCII file')
            return
        # n_s.k_version = '91bg'
        n_s.choose_model('EBV_model2', stype='st', RVhost=2.5)
        n_s.set_restbands()  # Auto pick appropriate rest-bands

        # Remove empty filters -- fix for 'ValueError: attempt to get argmin of an empty sequence'
        for filter in list(n_s.data.keys()):
            if len(n_s.data[filter].magnitude) == 0:
                del n_s.data[filter]
            elif self.origin == 'CSP' and filter in ['u', 'Y', 'J', 'H', 'Jrc2', 'Ydw']:
                print('[***] Special Process for CSP! Removing '+filter+'...')
                del n_s.data[filter]
        print('      Best filters:', list(n_s.data.keys()))

        for i in range(5):
            try:
                if self.origin == 'CSP':
                    initial_filters = []
                    for fil in ['B','V','g']:
                        if fil in list(n_s.data.keys()):
                            initial_filters.append(fil)
                    print('[***] Special Process for CSP! Fitting as '+str(initial_filters)+' -> remaining...')

                    n_s.fit(initial_filters, dokcorr=True, k_stretch=False, reset_kcorrs=True, **{'mangle': 1, 'calibration': 0})
                    n_s.fit(bands=None, dokcorr=True, k_stretch=False, reset_kcorrs=True, **{'mangle': 1, 'calibration': 0})
                else:
                    n_s.fit(bands=None, dokcorr=True, k_stretch=False, reset_kcorrs=True, **{'mangle': 1, 'calibration': 0})
                n_s.save(save_path)

                # Save paramaters
                for i in range(len(param_names)):
                    self.params.update({param_names[i]: {'value': n_s.parameters[snpy_param_names[i]],
                                                         'err': n_s.errors[snpy_param_names[i]]}})

                if show_plot:
                    n_s.plot(outfile=save_loc + 'plots/' + self.objname + '_snpyplots.png')
                    plt.show()
                    systime.sleep(3)
                plt.close()
                break
            except Exception as error:
                if 'All weights for filter' and 'are zero.' in str(error):
                    print('[!!!] Weights for filter', str(error).split(' ')[4], 'are zero. Removing...')
                    del n_s.data[str(error).split(' ')[4]]
                elif str(error) == 'Error:  to solve for EBVhost, you need to fit more than one filter':
                    print('[!!!] To few filters to fit!')
                    self.params.update({'mu': {'value': -1.0, 'err': -1.0}})
                    break
                else:
                    self.params.update({'mu': {'value': -1.0, 'err': -1.0}})
                    print(error)

        # Restore print statements
        sys.stdout = sys.__stdout__

        return
    def salt_fit(self, save_loc, show_plot=True, quiet=False):
        print('[+++] '+self.objname+' -- Fitting data with SALT3...')

        # Check quiet
        if quiet:
            sys.stdout = open(os.devnull, 'w')

        alpha, beta = float(CONSTANTS['salt_alpha']), float(CONSTANTS['salt_beta'])
        mB_const, M0 = float(CONSTANTS['salt_mB_const']), float(CONSTANTS['salt_absolute_mag'])
        try:
            # Fix filters
            filter_dict = {'u': 'cspu', 'g': 'cspg', 'r': 'cspr', 'i': 'cspi', 'B': 'cspB',
                           'V0': 'cspv3014', 'V1': 'cspv3009', 'V': 'cspv9844', 'Y': 'cspys',
                           'J': 'cspjs', 'Jrc2': 'cspjd', 'Jdw': 'cspjd', 'Ydw': 'cspyd', 'Hdw': 'csphd', 'H': 'csphs',
                           'c': 'atlasc', 'o': 'atlaso', 'ZTF_g': 'ztfg', 'ZTF_r': 'ztfr', 'ZTF_i': 'ztfi'}
            salt_time, salt_filters, salt_flux = np.array([]), np.array([]), np.array([])
            salt_dflux, salt_zp = np.array([]), np.array([])

            # Fix filters
            for i in range(len(self.filters)):
                if self.origin == 'CSP' and self.filters[i] in ['u', 'Y', 'J', 'H', 'Jrc2', 'Ydw']:
                    continue
                salt_time = np.append(salt_time, self.time[i])
                salt_filters = np.append(salt_filters, filter_dict[self.filters[i]])
                salt_flux = np.append(salt_flux, self.flux[i])
                salt_dflux = np.append(salt_dflux, self.dflux[i])
                salt_zp = np.append(salt_zp, self.zp[i])
            print('[***] Special Process for CSP!', np.unique(self.filters), '->', np.unique(salt_filters))

            data = Table([salt_time, salt_filters, salt_flux, salt_dflux, salt_zp, np.full(len(salt_time), 'ab')],
                         names=('time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'))

            # Create Model
            model = sncosmo.Model(source='salt3')

            # Fit data to model
            details = gen.TNS_details(self.coords[0], self.coords[1])
            model.set(z=self.z, t0=self.discovery_date)  # set the model's redshift.
            result, fitted_model = sncosmo.fit_lc(data, model, ['t0', 'x0', 'x1', 'c'], bounds={'x1': (-5, 5)})

            param_names = ['t0', 'x0', 'x1', 'c']
            for i in range(len(param_names)):
                self.params.update({param_names[i]: {'value': result.parameters[i+1],
                                                     'err': result.errors[param_names[i]]}})

            # Calculate
            pho_mB = -2.5 * np.log10(self.params['x0']['value']) + mB_const
            pho_mB_err = np.abs(-2.5 * (self.params['x0']['err'] / (self.params['x0']['value'] * np.log(10))))

            mu = pho_mB + (alpha * self.params['x1']['value']) - (beta * self.params['c']['value']) - M0
            mu_err = np.sqrt(pho_mB_err ** 2 + (np.abs(alpha) * self.params['x1']['err']) ** 2 + (np.abs(beta) * self.params['c']['err']) ** 2)

            self.params.update({'mu': {'value': mu, 'err': mu_err}})

            # Plot data with fit
            if show_plot:
                sncosmo.plot_lc(data, model=fitted_model, errors=result.errors)
                print('Saving plots to', save_loc + '/plots/' + self.objname + '_salt3lc.png')
                plt.savefig(save_loc + '/plots/' + self.objname + '_salt3lc.png')
                plt.show()

            print('      Successfully fit ' + self.objname + '!')
            print('Pausing for 1 seconds...')
            systime.sleep(1)
        except Exception as error:
            print(error)
            self.params.update({'mu': {'value': -1.0, 'err': -1.0}})

        # Restore print statements
        sys.stdout = sys.__stdout__

        return
    def get_host_mass(self, use_key=False):
        print('[+++] '+self.objname+' -- Finding host galaxy mass using GHOST...')
        local_coords = SkyCoord(self.coords[0], self.coords[1], unit=u.deg)
        galac_coords = local_coords.transform_to(Galactic())

        # Get CMB redshift
        helio_corr = (float(CONSTANTS['cmb_v_helio']) / float(CONSTANTS['cmb_c']) *
                      ((np.sin(galac_coords.b.deg) * np.sin(float(CONSTANTS['cmb_b_h'])) + np.cos(galac_coords.b.deg) *
                        np.cos(float(CONSTANTS['cmb_b_h'])) * np.cos(galac_coords.l.deg - float(CONSTANTS['cmb_l_h'])))))
        corr_term = 1 - helio_corr
        self.z_cmb = (1 + self.z) / corr_term - 1

        # Try mass key
        if use_key:
            mass_key = {}
            with open(CONSTANTS['mass_key_txt'], 'r') as f:
                temp = f.readlines()
                for line in temp:
                    line = line[:-1].split(', ')
                    if len(line) != 3:
                        continue
                    mass_key.update({line[0]: {'mass': line[1], 'mass_err': line[2]}})
            if self.objname in mass_key:
                print('      Found object in mass key! Pulling...')
                if float(mass_key[self.objname]['mass']) < 0:
                    print('      Found mass known to not be in GLADE/PanSTARR/SDSS')
                    return
                self.params.update({'hostMass': {'value': float(mass_key[self.objname]['mass']),
                                                 'err': float(mass_key[self.objname]['mass_err'])}})
                print('[+++] Mass taken from mass key!')
                return

        # Getting host data -- checks GLADE then PANSTARRS
        gMag, iMag, iAbsMaggMagErr, iMagErr, iAbsMagErr = -999.00, -999.00, -999.00, -999.00, -999.00
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                host_data = getTransientHosts(transientCoord=[local_coords],
                                              transientName=[self.objname],
                                              verbose=False,
                                              starcut="gentle", savepath='default/ghost_stuff/',
                                              GHOSTpath=CONSTANTS['ghost_data_loc'])
                print('      Identified Host Galaxy:', host_data.loc[0, 'NED_name'])

                # Get magnitudes from GLADE/PANSTARRS/SDSS
                if ~np.isnan(host_data.loc[0, 'gKronMag']) and ~np.isnan(host_data.loc[0, 'iKronMag']):
                    print('[+++] Mass taken from GLADE/PANSTARRS!')
                    gMag, iMag, iAbsMag = (host_data['gKronMag'].loc[0], host_data['iKronMag'].loc[0],
                                           host_data['iKronMag'].loc[0] - gen.current_cosmo().distmod(self.z).value)
                    gMagErr, iMagErr, iAbsMagErr = (host_data['gKronMagErr'].loc[0],
                                                    host_data['iKronMagErr'].loc[0], host_data['iKronMagErr'].loc[0])
        except Exception as error:
            print('Unknown GHOST error:', error)

        # Try SDSS
        if gMag == -999.00:
            print('[!!!] GHOST failed to find mass with GLADE/PANSTARRS, attempting to use SDSS...')
            result = SDSS.query_crossid(local_coords,
                                        photoobj_fields=['modelMag_g', 'modelMagErr_g', 'modelMag_i', 'modelMagErr_i'])
            if result == None:
                print('[!!!] GLADE/PANSTARRS/SDSS failed to find mass, returning zero mass...')
                self.params.update({'hostMass': {'value': 0.00, 'err': 0.00}})
                with open(CONSTANTS['mass_key_txt'], 'a') as f:
                    print('      Updating mass key with ' + self.objname + '...')
                    f.write(self.objname + ', -999.99, -999.99\n')
                return
            else:
                print('[+++] Mass taken from SDSS!')
                gMag, iMag, iAbsMag = (result['modelMag_g'].value[0], result['modelMag_i'].value[0],
                                       result['modelMag_i'].value[0] - gen.current_cosmo().distmod(self.z).value)
                gMagErr, iMagErr, iAbsMagErr = (result['modelMagErr_g'].value[0],
                                                result['modelMagErr_i'].value[0], result['modelMagErr_i'].value[0])

        # Mass Calculation -- Taylor et. al. 2011 -- eq. 8
        host_mass = (1.15 + (0.7 * (gMag - iMag)) - (0.4 * (iAbsMag)))

        # Error Propogation
        giMagErr = np.sqrt((gMagErr ** 2) + (iMagErr ** 2))
        host_mass_err = np.sqrt(((0.7 ** 2) * (giMagErr ** 2)) + ((0.4 ** 2) * (iAbsMagErr ** 2)))

        # Save Mass
        self.params.update({'hostMass': {'value': host_mass, 'err': host_mass_err}})
        print('      Success!', self.objname, 'host galaxy has a mass of:', host_mass, '+/-', host_mass_err, 'log(M_*/[M_sun])')

        # Update mass key
        if use_key:
            with open(CONSTANTS['mass_key_txt'], 'a') as f:
                print('      Updating mass key with ' + self.objname + '...')
                f.write(self.objname + ', ' + str(host_mass) + ', ' + str(host_mass_err) + '\n')

        print('      Removing GHOST data...')
        shutil.rmtree('default/ghost_stuff/')  # Clear messy data
        os.mkdir('default/ghost_stuff/')
        return
# ------------------------------------------------------------------------------------------------------------------- #
def combined_fit(algo='snpy', dmag_max=0, dflux_max=0):
    sys.stdout = open(os.devnull,'w') # Lots of unecessary output

    # Load and set up files
    csp_files, atlas_files, ztf_files = glob.glob('data/CSP/*.txt'), glob.glob('data/ATLAS/*.txt'), glob.glob('data/ZTF/*.txt')
    CSP_SNe, ATLAS_SNe, ZTF_SNe, atlas_names = {}, {}, {}, []
    for file in csp_files:
        tempSN = sn91bg().make_class('CSP', file, dmag_max, dflux_max)
        if tempSN is not None:
            CSP_SNe.update({tempSN.objname: tempSN})
    for file in atlas_files:
        tempSN = sn91bg().make_class('ATLAS', file, dmag_max, dflux_max)
        if tempSN is not None:
            ATLAS_SNe.update({tempSN.objname: tempSN})
    for file in ztf_files:
        tempSN = sn91bg().make_class('ZTF', file, dmag_max, dflux_max)
        if tempSN is not None:
            ZTF_SNe.update({tempSN.objname: tempSN})

    sys.stdout = sys.__stdout__

    # List of overlapp
    altas_ztf_list = []
    for element in list(ATLAS_SNe.keys()):
        if element in list(ZTF_SNe.keys()):
            altas_ztf_list.append(element)

    # List of all uniuqe SNe
    combined_list = np.unique(np.hstack((np.hstack((list(CSP_SNe.keys()), list(ATLAS_SNe.keys()))),
                                         list(ZTF_SNe.keys()))))

    # Combine filters and data
    combined_SNe = []
    for name in combined_list:
        if name in altas_ztf_list:
            new_zp = np.hstack((ATLAS_SNe[name].zp, ZTF_SNe[name].zp))
            new_filters = np.hstack((ATLAS_SNe[name].filters, ZTF_SNe[name].filters))
            new_time = np.hstack((ATLAS_SNe[name].time, ZTF_SNe[name].time))
            new_flux = np.hstack((ATLAS_SNe[name].flux, ZTF_SNe[name].flux))
            new_dflux = np.hstack((ATLAS_SNe[name].dflux, ZTF_SNe[name].dflux))
            new_mag = np.hstack((ATLAS_SNe[name].mag, ZTF_SNe[name].mag))
            new_dmag = np.hstack((ATLAS_SNe[name].dmag, ZTF_SNe[name].dmag))

            ATLAS_SNe[name].origin = 'ATLAS-ZTF'
            ATLAS_SNe[name].zp = new_zp
            ATLAS_SNe[name].filters = new_filters
            ATLAS_SNe[name].time = new_time
            ATLAS_SNe[name].flux = new_flux
            ATLAS_SNe[name].dflux = new_dflux
            ATLAS_SNe[name].mag = new_mag
            ATLAS_SNe[name].dmag = new_dmag
            combined_SNe.append(ATLAS_SNe[name])
        elif name in list(CSP_SNe.keys()):
            combined_SNe.append(CSP_SNe[name])
        elif name in list(ATLAS_SNe.keys()):
            combined_SNe.append(ATLAS_SNe[name])
        elif name in list(ZTF_SNe.keys()):
            combined_SNe.append(ZTF_SNe[name])

    # Fitting data
    fit_combined_SNe = []
    for SN in combined_SNe:
        print('-----------------------------------------------------------------------------------------------')
        print('[', combined_SNe.index(SN)+1, '/', len(combined_SNe), '] Fitting data for '+SN.objname+' ['+algo+'] ['+SN.origin+']...')
        SN.fit(algo=algo, save_loc=CONSTANTS[algo+'_combined_saved_loc'])
        # Check if fit failed
        if (SN != None and
            SN.params['mu']['value'] > 0 and
            'hostMass' in SN.params and
            SN.params['hostMass']['value'] > 0):
            fit_combined_SNe.append(SN)
    print('Sucessfully fit [', len(fit_combined_SNe), '/', len(combined_SNe), ']!')

    # Save parameters to file
    save_params_to_file(CONSTANTS[algo+'_combined_saved_loc']+'combined_params.txt', fit_combined_SNe)

    # Cut sample
    sample_cutter(CONSTANTS[algo+'_combined_saved_loc']+'combined_params.txt', algo=algo)

    return fit_combined_SNe
def batch_fit(data_set, algo='snpy', dmag_max=0, dflux_max=0):
    SNe, files = [], glob.glob(CONSTANTS[data_set.lower()+'_data_loc'] + '*.txt')
    for path in files:
        print('[', files.index(path) + 1, '/', len(files), ']')
        print('-----------------------------------------------------------------------------------------------')
        tempSN = sn91bg().make_class(data_set=data_set, path=path, dmag_max=dmag_max, dflux_max=dflux_max)
        tempSN.fit(algo)

        # Check if fit failed
        if (tempSN != None and
            tempSN.params['mu']['value'] > 0 and
            'hostMass' in tempSN.params and
            tempSN.params['hostMass']['value'] > 0):
            SNe.append(tempSN)
    print('Sucessfully fit [', len(SNe), '/', len(files), ']!')

    # Save params to file
    save_params_to_file(CONSTANTS[algo+'_'+data_set.lower() + '_saved_loc'] + data_set.lower() + '_params.txt', SNe)

    # Cut sample
    sample_cutter(CONSTANTS[algo+'_'+data_set.lower() + '_saved_loc'] + data_set.lower() + '_params.txt', algo=algo)

    return SNe
def sample_cutter(path, algo='snpy'):
    data = np.genfromtxt(path, delimiter=', ', skip_header=1, dtype=str)
    original_num = str(len(data[:, 0]))
    if len(data) == 0:
        print('[+++] No params to cut!')
        return
    with open(path, 'r') as f:
        hdr = f.readline().split(', ')
        hdr[-1] = hdr[-1][:-1]
    if algo == 'snpy':
        print('[+++] Cutting sample for SNooPy data...')
        # cuts = {'z': 0.015, 'EBVhost': (-1, 1), 'EBVhost_err': 0.2, 'st': (0.3, 1.0), 'st_err': 0.1, 'Tmax_err': 1} # Old
        cuts = {'z': 0.015, 'EBVhost': (-1, 0.28), 'EBVhost_err': 0.1, 'st': (0.3, 1.0), 'st_err': 0.1, 'Tmax_err': 1} # Maximized
        f_out = '      | '
        for c in cuts:
            f_out += c + ': ' + str(cuts[c]) + ' | '
        print(f_out)
        data = data[(data[:, hdr.index('z_cmb')].astype(float) > cuts['z']) &
                    (data[:, hdr.index('EBVhost')].astype(float) > cuts['EBVhost'][0]) &
                    (data[:, hdr.index('EBVhost')].astype(float) < cuts['EBVhost'][1]) &
                    (data[:, hdr.index('EBVhost_err')].astype(float) < cuts['EBVhost_err']) &
                    (data[:, hdr.index('st')].astype(float) > cuts['st'][0]) &
                    (data[:, hdr.index('st')].astype(float) < cuts['st'][1]) &
                    (data[:, hdr.index('st_err')].astype(float) < cuts['st_err']) &
                    (data[:, hdr.index('Tmax_err')].astype(float) < cuts['Tmax_err'])]
    elif algo == 'salt':
        print('[+++] Cutting sample for SALT data...')
        # cuts = {'z': 0.015, 'x0': 999, 'x0_err': 999, 'x1': (-5, 5), 'x1_err': 1, 'c': (-0.3, 999), 'c_err': 0.1, 't0_err': 2}
        cuts = {'z': 0.015, 'x0': 999, 'x0_err': 999, 'x1': (-4.5, 4.5), 'x1_err': 1, 'c': (-0.3, 1), 'c_err': 0.3, 't0_err': 0.5}
        f_out = '      | '
        for c in cuts:
            f_out += c + ': ' + str(cuts[c]) + ' | '
        print(f_out)
        data = data[(data[:, hdr.index('z_cmb')].astype(float) > cuts['z']) &
                    (data[:, hdr.index('x0')].astype(float) < cuts['x0']) &
                    (data[:, hdr.index('x0_err')].astype(float) < cuts['x0_err']) &
                    (data[:, hdr.index('x1')].astype(float) > cuts['x1'][0]) &
                    (data[:, hdr.index('x1')].astype(float) < cuts['x1'][1]) &
                    (data[:, hdr.index('x1_err')].astype(float) < cuts['x1_err']) &
                    (data[:, hdr.index('c')].astype(float) > cuts['c'][0]) &
                    (data[:, hdr.index('c')].astype(float) < cuts['c'][1]) &
                    (data[:, hdr.index('c_err')].astype(float) < cuts['c_err']) &
                    (data[:, hdr.index('t0_err')].astype(float) < cuts['t0_err'])]

    # Save to file
    with open(path[:-4]+'_cut.txt', 'w') as f:
        f_out = hdr[0]
        for h in hdr[1:]:
            f_out += ', ' + h
        f.write(f_out + '\n')
        for line in data[:]:
            f_out = line[0]
            for item in line[1:]:
                f_out +=  ', ' + str(item)
            f.write(f_out + '\n')
    print('      Cut file saved to...', path[:-4]+'_cut.txt')
    print('      [ '+str(len(data[:, 0]))+' / '+original_num+' ]')
    return
# ------------------------------------------------------------------------------------------------------------------- #
def batch_load(data_set, algo='snpy'):
    SNe = []
    for path in glob.glob('saved/' + algo + '/' + data_set.lower() + '/classes/*_class.txt'):
        SNe.append(sn91bg().load_from_file(path))
    return SNe
# ------------------------------------------------------------------------------------------------------------------- #
def save_params_to_file(save_loc, SNe):
    print('[+++] Saving params to '+save_loc+'...')
    if len(SNe) == 0:
        print('[+++] No params to save!')
        return
    with open(save_loc, 'w') as f:
        hdr = 'objname, ra, dec, z, z_cmb, MJDs, MJDe, origin'
        for params in SNe[0].params:
            hdr += ', ' + str(params) + ', ' + str(params) + '_err'
        f.write(hdr + '\n')
        for SN in SNe:
            line = (SN.objname + ', ' + str(SN.coords[0]) + ', ' + str(SN.coords[1]) + ', ' +
                    str(SN.z) + ', ' + str(SN.z_cmb) + ', ' +
                    str(SN.period[0]) + ', ' + str(SN.period[1]) + ', ' + str(SN.origin))
            for param in SN.params:
                line += ', ' + str(SN.params[param]['value']) + ', ' + str(SN.params[param]['err'])
            f.write(line + '\n')
    return
# ------------------------------------------------------------------------------------------------------------------- #
def residual_plotter(path, x_params='z', data_set='', algo='', labels=False, save_loc=None):
    # Pull data from saved text
    data = np.genfromtxt(path, delimiter=', ', skip_header=1, dtype=str)
    if len(data_set) == 0:
        data_set = path.split('/')[-1].split('_')[0].upper()
    if len(algo) == 0:
        algo = path.split('/')[-3].upper()
    if save_loc is None:
        save_loc = path[:-len(path.split('/')[-1])]+data_set.lower()+'_resid_v_'+x_params.lower()+'.png'

    # Get header
    with open(path, 'r') as f:
        hdr = f.readline().split(', ')
        hdr[-1] = hdr[-1][:-1]

    fig, axs = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'width_ratios': [10, 1]}, constrained_layout=True)
    color_wheel = {'ZTF': '#81ADC8', 'ATLAS': '#EEDFAA', 'CSP': '#CD4631', 'ATLAS-ZTF': '#DEA47E',
                   'ZTF_SNPY': '#FFADC8', 'ATLAS_SNPY': '#FFDFAA', 'CSP_SNPY': '#FF4631', 'ATLAS-ZTF_SNPY': '#FFA47E',
                   'ZTF_SALT': '#AAADC8', 'ATLAS_SALT': '#AADFAA', 'CSP_SALT': '#AA4631', 'ATLAS-ZTF_SALT': '#AAA47E'}

    # Determine x-axis
    if x_params == 'mass':
        axs[0].set(xlabel=' Host Stellar Mass (log $M_{*}$/$[M_{\odot}]$)', ylabel='Hubble Residuals (mag)')  # Sub-plot Labels
        # axs[0].set_xlim(7.5, 12); axs[0].set_ylim(-3.0, 3.0); axs[1].set_ylim(-3.0, 3.0)
        x_axis, x_axis_err = data[:, hdr.index('hostMass')].astype(float), data[:, hdr.index('hostMass_err')].astype(float)
        title = "Hubble Residuals vs. Host Stellar Mass of '"+data_set+"' 91bg-like SNe Ia\n"
    elif x_params == 'z':
        axs[0].set(xlabel='Host Galaxy CMB Redshift', ylabel='Hubble Residuals (mag)')  # Sub-plot Labels
        # axs[0].set_xlim(0.01, 0.08); axs[0].set_ylim(-4.0, 4.0); axs[1].set_ylim(-4.0, 4.0)
        x_axis, x_axis_err = data[:, hdr.index('z_cmb')].astype(float), np.full(len(data[:, hdr.index('z_cmb')]), np.nan)
        title = "Hubble Residuals vs. Redshift of '"+data_set+"' 91bg-like SNe Ia\n"
    else:
        raise ValueError("[!!!] Invalid x_params ['mass'/'z']")

    # Plot points
    for origin in np.unique(data[:, hdr.index('origin')]):
        indexs = np.where(data[:, hdr.index('origin')] == origin)[0]
        resid_mu = data[:, hdr.index('mu')].astype(float)[indexs] - gen.current_cosmo().distmod(data[:, hdr.index('z_cmb')].astype(float)[indexs]).value
        resid_mu_err = data[:, hdr.index('mu_err')].astype(float)[indexs]

        x_axis_plt, x_axis_err_plt = x_axis[indexs], x_axis_err[indexs]

        axs[0].errorbar(x_axis_plt, resid_mu, xerr=x_axis_err_plt, yerr=resid_mu_err,
                        color=color_wheel[origin], fmt='o', label=origin)

        # Labels
        if labels:
            for i in range(len(resid_mu)):
                axs[0].text(x_axis_plt[i], resid_mu[i], str(data[:, 0][indexs][i]), size='x-small', va='top')

    # Plot histogram
    all_resid_mu = data[:, hdr.index('mu')].astype(float) - gen.current_cosmo().distmod(data[:, hdr.index('z_cmb')].astype(float)).value
    all_resid_mu_err = data[:, hdr.index('mu_err')].astype(float)
    axs[1].hist(all_resid_mu, bins=15, orientation="horizontal", color='#9E6240')

    # Adjusting axies presentation
    ylimiter = np.max(np.abs(all_resid_mu))+0.3
    axs[0].set_ylim(-ylimiter, ylimiter); axs[1].set_ylim(-ylimiter, ylimiter)

    # Formatting
    fig.suptitle(title +
                 'Scatter: ' + str(round(np.std(all_resid_mu), 2)) + ' | # of pts: ' + str(len(all_resid_mu)) +
                 ' | SNR: ' + str(round(np.sqrt(np.abs(np.average(all_resid_mu)) / np.abs(np.std(all_resid_mu_err))), 2)),
                 size='medium')
    axs[1].get_yaxis().set_visible(False) # Turn off y-axis labels
    axs[0].legend(loc='best')
    print('Saved figure to... ', save_loc)
    plt.savefig(save_loc)
    plt.show()
    return
def histogram_plotter(path, param_bins=[45, 45, 45, 45, 45]):
    # Pull data from saved text
    data = np.genfromtxt(path, delimiter=', ', skip_header=1, dtype=str)
    data_set = path.split('/')[-1].split('_')[0].upper()
    algo = path.split('/')[-3].upper()

    # Pull data
    if algo == 'SNPY':
        mu, st, Tmax, EBVhost, hostMass = (data[:, 8].astype(float), data[:, 10].astype(float), data[:, 12].astype(float),
                                            data[:, 14].astype(float), data[:, 16].astype(float))
        params = [mu, st, Tmax, EBVhost, hostMass]
    elif algo == 'SALT':
        t0, x0, x1, c, mu, hostMass = (data[:, 8].astype(float), data[:, 10].astype(float), data[:, 12].astype(float),
                                       data[:, 14].astype(float), data[:, 16].astype(float), data[:, 18].astype(float))
        params = [t0, x0, x1, c, mu, hostMass]

    # Plot
    fig, ax = plt.subplots(1, 5, figsize=(16, 4), layout='constrained')
    param_names = ['mu', 'st', 'Tmax', 'EBVhost', 'host_mass']
    for i in range(len(params)):
        ax[i].hist(params[i], bins=param_bins[i])
        if i != 0:
            ax[i].get_yaxis().set_visible(False)
        ax[i].set_xlabel(param_names[i])

    plt.suptitle("Parameters for '" + path.split('/')[-1].split('_')[0].upper()
                 + "' data\n Number of Transients: " + str(len(data)), fontsize=20)
    print('Saved figure to... ', path.split('/')[-1].split('_')[0]+'_hist.png')
    plt.savefig(path.split('/')[-1].split('_')[0]+'_hist.png')
    plt.show()
    return
def make_display_plots():
    residual_plotter('output/COMBINED_combined_snpy_params_cut.txt',
                     x_params='z', data_set='CSP-ATLAS-ZTF (snpy)', algo='snpy',
                     save_loc='saved/readme_plots/csp-atlas-ztf_snpy_resid_v_z.png')
    residual_plotter('output/COMBINED_combined_salt_params_cut.txt',
                     x_params='z', data_set='CSP-ATLAS-ZTF (salt)', algo='salt',
                     save_loc='saved/readme_plots/csp-atlas-ztf_salt_resid_v_z.png')
    residual_plotter('output/merged_params_cut.txt',
                     x_params='z', data_set='CSP-ATLAS-ZTF (snpy-salt)', algo='snpy-salt',
                     save_loc='saved/readme_plots/merged_resid_v_z.png')

    mass_step_residual_plot('output/COMBINED_combined_snpy_params_cut.txt',
                            data_set='CSP-ATLAS-ZTF (snpy)', algo='snpy',
                            save_loc='saved/readme_plots/csp-atlas-ztf_snpy_resid_v_mass.png')
    mass_step_residual_plot('output/COMBINED_combined_salt_params_cut.txt',
                            data_set='CSP-ATLAS-ZTF (salt)', algo='salt',
                            save_loc='saved/readme_plots/csp-atlas-ztf_salt_resid_v_mass.png')
    mass_step_residual_plot('output/merged_params_cut.txt',
                            data_set='CSP-ATLAS-ZTF (snpy-salt)', algo='snpy-salt',
                            save_loc='saved/readme_plots/merged_resid_v_mass.png')

    return
def param_hist_compare(set1, set2, bin_width, xrange=None, title=''):
    # Select bins
    if bin_width == 'auto':
        bin1, bin2 = None, None
    elif len(bin_width) > 1:
        bin1, bin2 = int((np.max(set1) - np.min(set1)) / bin_width[0]), int((np.max(set2) - np.min(set2)) / bin_width[1])
    else:
        bin1, bin2 = int((np.max(set1) - np.min(set1)) / bin_width), int((np.max(set2) - np.min(set2)) / bin_width)

    # Plotting
    plt.figure(figsize=(12, 6))
    plt.hist(set1, color='#62BEC1', label=title.split(' ')[0], alpha=1,
             bins=bin1)
    plt.hist(set2, color='#5AD2F4', label=title.split(' ')[2], alpha=0.825,
             bins=bin2)

    # Data Lines
    plt.axvline(x=np.median(set1), label='Median '+title.split(' ')[0], linewidth=4, color='#52a1a3', linestyle='--')
    plt.axvline(x=np.median(set2), label='Median '+title.split(' ')[2], linewidth=4, color='#4bb0cc', linestyle='-.')

    # Formatting
    plt.title(title)
    plt.xlim(xrange)
    plt.legend()

    # Save Plot
    plt.savefig('saved/HiCATvPan+_'+title.split(' ')[4].lower()+'.png')
    print('Saved to...', 'saved/HiCATvPan+_'+title.split(' ')[4].lower()+'.png')

    plt.show()
    return
def mass_step_residual_plot(path, data_set='', title='', labels=False, save_loc=None):
    # Pull data from saved text
    data = np.genfromtxt(path, delimiter=', ', skip_header=1, dtype=str)
    if len(data_set) == 0:
        data_set = 'default'
    if save_loc is None:
        save_loc = path[:-len(path.split('/')[-1])]+data_set.lower()+'_resid_v_mass.png'
    if len(title) == 0:
        title = "Hubble Residuals vs. Host Stellar Mass of '" + data_set + "' 91bg-like SNe Ia\n"

    # Get header
    with open(path, 'r') as f:
        hdr = f.readline().split(', ')
        hdr[-1] = hdr[-1][:-1]

    fig, axs = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'width_ratios': [10, 1]}, constrained_layout=True)
    color_wheel = {'ZTF': '#81ADC8', 'ATLAS': '#EEDFAA', 'CSP': '#CD4631', 'ATLAS-ZTF': '#DEA47E',
                   'ZTF_SNPY': '#FFADC8', 'ATLAS_SNPY': '#FFDFAA', 'CSP_SNPY': '#FF4631', 'ATLAS-ZTF_SNPY': '#FFA47E',
                   'ZTF_SALT': '#AAADC8', 'ATLAS_SALT': '#AADFAA', 'CSP_SALT': '#AA4631', 'ATLAS-ZTF_SALT': '#AAA47E',
                   'Pan+': 'maroon'}

    # Determine x-axis
    axs[0].set(xlabel=' Host Stellar Mass (log $M_{*}$/$[M_{\odot}]$)', ylabel='Hubble Residuals (mag)')  # Sub-plot Labels
    x_axis, x_axis_err = data[:, hdr.index('hostMass')].astype(float), data[:, hdr.index('hostMass_err')].astype(float)

    # Plot points
    all_mu, all_mu_err = data[:, hdr.index('mu')].astype(float), data[:, hdr.index('mu_err')].astype(float)
    all_z = data[:, hdr.index('z_cmb')].astype(float)
    all_mass = data[:, hdr.index('hostMass')].astype(float)

    for origin in np.unique(data[:, hdr.index('origin')]):
        indexs = np.where(data[:, hdr.index('origin')] == origin)[0]
        resid_mu = data[:, hdr.index('mu')].astype(float)[indexs] - gen.current_cosmo().distmod(data[:, hdr.index('z_cmb')].astype(float)[indexs]).value
        resid_mu_err = data[:, hdr.index('mu_err')].astype(float)[indexs]

        x_axis_plt, x_axis_err_plt = x_axis[indexs], 0.00

        axs[0].errorbar(x_axis_plt, resid_mu, xerr=x_axis_err_plt, yerr=resid_mu_err,
                        color=color_wheel[origin], fmt='o', label=origin, alpha=1)

        # Labels
        if labels:
            for i in range(len(resid_mu)):
                axs[0].text(x_axis_plt[i], resid_mu[i], str(data[:, 0][indexs][i]), size='x-small', va='top')

    # Plot histogram
    all_resid_mu = all_mu - gen.current_cosmo().distmod(all_z).value
    all_resid_mu_err = np.copy(all_mu_err)
    axs[1].hist(all_resid_mu, bins=30, orientation="horizontal", color='#9E6240')

    # Plot mass step lines
    cut = round(np.median(all_mass), 4)
    lower_resid = np.average(all_mu[all_mass < cut] - gen.current_cosmo().distmod(all_z[all_mass < cut]).value,
                             weights=1/(all_mu_err[all_mass < cut]**2))
    upper_resid = np.average(all_mu[all_mass > cut] - gen.current_cosmo().distmod(all_z[all_mass > cut]).value,
                             weights=1/(all_mu_err[all_mass > cut]**2))
    lower_resid_err = np.std(all_mu_err[all_mass < cut]) / np.sqrt(len(all_mu_err[all_mass < cut]))
    upper_resid_err = np.std(all_mu_err[all_mass > cut]) / np.sqrt(len(all_mu_err[all_mass > cut]))
    mass_step = np.abs(lower_resid - upper_resid)
    mass_step_err = np.sqrt((lower_resid_err**2) + (upper_resid_err**2))

    # Plot before step
    axs[0].hlines(y=lower_resid, xmin=np.min(all_mass)-0.1, xmax=cut, linestyle='--', linewidth=2.5, color='k')
    axs[0].fill_between([np.min(all_mass)-0.1, cut], lower_resid - lower_resid_err, lower_resid + lower_resid_err, color='k', alpha=0.2)

    # Plot after step
    axs[0].hlines(y=upper_resid, xmin=cut, xmax=np.max(all_mass)+0.1, linestyle='--', linewidth=2.5, color='k')
    axs[0].fill_between([cut, np.max(all_mass)+0.1], upper_resid - upper_resid_err, upper_resid + upper_resid_err, color='k', alpha=0.2)

    # Plot vertical line between points
    axs[0].vlines(x=cut, ymin=lower_resid, ymax=upper_resid, linestyle='--', linewidth=2.5, color='k')

    # Adjusting axies presentation
    ylimiter = np.max(np.abs(all_resid_mu))+0.3
    axs[0].set_ylim(-ylimiter, ylimiter); axs[1].set_ylim(-ylimiter, ylimiter)

    # Formatting
    fig.suptitle(title +
                 'Scatter: ' + str(round(np.std(all_resid_mu), 2)) + ' | # of pts: ' + str(len(all_resid_mu)) +
                 ' | SNR: ' + str(round(np.sqrt(np.abs(np.average(all_resid_mu)) / np.abs(np.std(all_resid_mu_err))), 2)) +
                 '\nMass Step ('+str(cut)+'): ' + str(round(mass_step, 4)) + ' +/- ' + str(round(mass_step_err, 6)),
                 size='medium')
    axs[1].get_yaxis().set_visible(False) # Turn off y-axis labels
    axs[0].legend(loc='best')
    print('Saved figure to... ', save_loc)
    plt.savefig(save_loc)
    plt.show()
    return
# ------------------------------------------------------------------------------------------------------------------- #
def mass_step_calc(path, cut=10, quiet=False):
    data = np.genfromtxt(path, delimiter=', ', skip_header=1, dtype=str)
    with open(path, 'r') as f:
        hdr = f.readline().split(', ')
        hdr[-1] = hdr[-1][:-1]

    source = path.split('/')[-1].split('.')[0].upper().replace('_', ' ')+' '+str(cut).upper()
    mu, mu_err = data[:, hdr.index('mu')].astype(float), data[:, hdr.index('mu_err')].astype(float)
    mass = data[:, hdr.index('hostMass')].astype(float)
    z = data[:, hdr.index('z_cmb')].astype(float)

    if cut == 'median':
        cut = round(np.median(mass), 4)

    all_resid = mu - gen.current_cosmo().distmod(z).value
    lower_resid, lower_resid_err = mu[mass < cut] - gen.current_cosmo().distmod(z[mass < cut]).value, mu_err[mass < cut]
    upper_resid, upper_resid_err = mu[mass > cut] - gen.current_cosmo().distmod(z[mass > cut]).value, mu_err[mass > cut]

    mass_step = np.abs(np.average(lower_resid) - np.average(upper_resid))
    mass_step_err = np.abs(np.average(lower_resid_err) - np.average(upper_resid_err))

    if quiet:
        sys.stdout = open(os.devnull, 'w')

    print('***********************************************************************************************************')
    print(source)
    print('Mass Step [M(<' + str(cut) + ') - M(>' + str(cut) + ')]:', round(mass_step, 4), '+/-', round(mass_step_err, 4))
    print('\t----------------------')
    print('\tScatter:', round(np.std(all_resid), 4),
          '| M(<' + str(cut) + '):', round(np.std(lower_resid), 4),
          '| M(>' + str(cut) + '):', round(np.std(upper_resid), 4))
    print('\t----------------------')
    print('\tNumber of Targets:', len(mu),
          '| M(<' + str(cut) + '):', len(mu[mass < cut]),
          '| M(>' + str(cut) + '):', len(mu[mass > cut]))

    sys.stdout = sys.__stdout__

    return mass_step, mass_step_err, {'cut': cut, 'scatter': np.std(all_resid), 'source': source}
def alt_mass_step_calc(mu, mu_err, mass, z, cut=10, quiet=False):
    if cut == 'median':
        cut = round(np.median(mass), 4)

    all_resid = mu - gen.current_cosmo().distmod(z).value

    lower_resid = np.average(mu[mass < cut] - gen.current_cosmo().distmod(z[mass < cut]).value,
                             weights=1/(mu_err[mass < cut]**2))
    upper_resid = np.average(all_mu[mass > cut] - gen.current_cosmo().distmod(z[mass > cut]).value,
                             weights=1/(mu_err[mass > cut]**2))
    lower_resid_err = np.std(mu_err[mass < cut]) / np.sqrt(len(mu_err[mass < cut]))
    upper_resid_err = np.std(mu_err[mass > cut]) / np.sqrt(len(mu_err[mass > cut]))
    mass_step = np.abs(lower_resid - upper_resid)
    mass_step_err = np.sqrt((lower_resid_err**2) + (upper_resid_err**2))

    if quiet:
        sys.stdout = open(os.devnull, 'w')

    print('***********************************************************************************************************')
    print('Mass Step [M(<' + str(cut) + ') - M(>' + str(cut) + ')]:',
          round(mass_step, 4), '+/-', round(mass_step_err, 4),
          '[', round(mass_step-mass_step_err, 4), '-', round(mass_step+mass_step_err, 4), ']')
    print('\t----------------------')
    print('\tScatter:', round(np.std(all_resid), 4),
          '| M(<' + str(cut) + '):', round(np.std(lower_resid), 4),
          '| M(>' + str(cut) + '):', round(np.std(upper_resid), 4))
    print('\t----------------------')
    print('\tNumber of Targets:', len(mu),
          '| M(<' + str(cut) + '):', len(mu[mass < cut]),
          '| M(>' + str(cut) + '):', len(mu[mass > cut]))

    sys.stdout = sys.__stdout__

    return mass_step, mass_step_err, {'cut': cut, 'scatter': np.std(all_resid)}
def merge_snpy_salt(snpy_path, salt_path, save_loc):
    snpy_data = np.genfromtxt(snpy_path, delimiter=', ', skip_header=1, dtype=str)
    salt_data = np.genfromtxt(salt_path, delimiter=', ', skip_header=1, dtype=str)

    with open(snpy_path, 'r') as f:
        hdr_snpy = f.readline().split(', ')
        hdr_snpy[-1] = hdr_snpy[-1][:-1]
    with open(salt_path, 'r') as f:
        hdr_salt = f.readline().split(', ')
        hdr_salt[-1] = hdr_salt[-1][:-1]

    source, n = ['snpy', 'salt'], 0
    names = []
    with open(save_loc, 'w') as f:
        print('objname, z_cmb, origin, mu, mu_err, hostMass, hostMass_err', file=f)
        for data in [snpy_data, salt_data]:
            for i in range(len(data[:, 0])):
                if data[i, 0] not in names:
                    if source[n] == 'snpy':
                        print(data[i, hdr_snpy.index('objname')] + ',',
                              data[i, hdr_snpy.index('z_cmb')] + ',',
                              data[i, hdr_snpy.index('origin')] + '_' + source[n].upper() + ',',
                              data[i, hdr_snpy.index('mu')] + ',',
                              data[i, hdr_snpy.index('mu_err')] + ',',
                              data[i, hdr_snpy.index('hostMass')] + ',',
                              data[i, hdr_snpy.index('hostMass_err')],
                              file=f)
                    elif source[n] == 'salt':
                        print(data[i, hdr_salt.index('objname')] + ',',
                              data[i, hdr_salt.index('z_cmb')] + ',',
                              data[i, hdr_salt.index('origin')] + '_' + source[n].upper() + ',',
                              data[i, hdr_salt.index('mu')] + ',',
                              data[i, hdr_salt.index('mu_err')] + ',',
                              data[i, hdr_salt.index('hostMass')] + ',',
                              data[i, hdr_salt.index('hostMass_err')],
                              file=f)
                    names.append(data[i, 0])
            n += 1
    print('Final number of objects:', len(names),
          '[snpy: '+str(len(snpy_data[:, 0]))+']',
          '[salt: '+str(len(salt_data[:, 0]))+']')
    print('Saved merged param file to... ', save_loc)

    return
def mass_step_scatter_comparison():
    mass_steps, mass_step_errs, scatters = [], [], []
    print('source, mass_step, mass_step_err, params')
    paths = ['saved/salt/combined/combined_params_cut.txt',
             'saved/salt/combined/combined_params.txt',
             'saved/snpy/combined/combined_params_cut.txt',
             'saved/snpy/combined/combined_params.txt',
             'saved/snpy_salt_combined_params_cut.txt',
             'saved/snpy_salt_combined_params.txt']
    labels = ['salt_cut_10',
              'salt_cut_median',
              'salt_uncut_10',
              'salt_uncut_median',
              'snpy_cut_10',
              'snpy_cut_median',
              'snpy_uncut_10',
              'snpy_uncut_median',
              'merged_cut_10',
              'merged_cut_median',
              'merged_uncut_10',
              'merged_uncut_median']
    cut = 'median'
    for path in paths:
        for cut in [10, 'median']:
            mass_step, mass_step_err, params = mass_step_calc(path, cut=cut, quiet=True)
            # print(path.split('/')[-1].split('.')[0].upper().replace('_', ' ')+', ',
            #       str(round(mass_step, 7))+', ', str(round(mass_step_err, 7))+', ', params)
            mass_steps.append(mass_step)
            mass_step_errs.append(mass_step_err)
            scatters.append(params['scatter'])

    plt.figure(figsize=(12, 6))
    plt.title('Mass Step v. Scatter \nCut = '+str(cut)+' log $M_{*}$/$[M_{\odot}]$')
    for i in range(len(mass_steps)):
        plt.errorbar(mass_steps[i], scatters[i], xerr=mass_step_errs[i], fmt='o')
        plt.text(mass_steps[i], scatters[i], labels[i], size='x-small', va='top')
    plt.ylabel('Scatter'); plt.xlabel('Mass Step')
    plt.show()

    with open('mass_step_scatter.txt', 'w') as f:
        f.write('source, mass_step, mass_step_err, scatter\n')
        for i in range(len(mass_steps)):
            f.write(str(labels[i])+', '+str(mass_steps[i])+', '+str(mass_step_errs[i])+', '+str(scatters[i])+'\n')

    return
def mass_step_compare_to_normals():
    hicat_salt_data = np.genfromtxt('output/COMBINED_combined_salt_params_cut.txt', delimiter=', ', skip_header=1,
                                    dtype=str)
    hicat_snpy_data = np.genfromtxt('output/COMBINED_combined_snpy_params_cut.txt', delimiter=', ', skip_header=1,
                                    dtype=str)
    hicat_merged_data = np.genfromtxt('saved/snpy_salt_combined_params_cut.txt', delimiter=', ', skip_header=1,
                                      dtype=str)
    panplus_data = np.genfromtxt('txts/PanPlus_Latest.FITRES', skip_header=69, dtype=str)[:,
                   1:]  # Removes 'SN:' at begining of each line
    hdr = ('CID CIDint IDSURVEY TYPE FIELD CUTFLAG_SNANA zHEL zHELERR zCMB zCMBERR zHD zHDERR VPEC VPECERR MWEBV '
           'HOST_LOGMASS HOST_LOGMASS_ERR HOST_sSFR HOST_sSFR_ERR PKMJDINI SNRMAX1 SNRMAX2 SNRMAX3 PKMJD PKMJDERR '
           'x1 x1ERR c cERR mB mBERR x0 x0ERR COV_x1_c COV_x1_x0 COV_c_x0 NDOF FITCHI2 FITPROB RA DEC HOST_RA HOST_DEC '
           'HOST_ANGSEP TGAPMAX TrestMIN TrestMAX ELU HOSTGAL_SFR HOSTGAL_SFR_ERR HOSTGAL_sSFR HOSTGAL_sSFR_ERR CUTMASK '
           'MU MUMODEL MUERR MUERR_RENORM MUERR_RAW MUERR_VPEC MURES MUPULL M0DIF M0DIFERR CHI2 biasCor_mu biasCorErr_mu '
           'biasCor_mB biasCor_x1 biasCor_c biasScale_muCOV IDSAMPLE IZBIN').split(' ')
    panplus_data = panplus_data[panplus_data[:, hdr.index('zCMB')].astype(float) < 0.08]  # Only low-z SNe

    print('\nMass Step of Norm Ia Panth+')
    norm_step, norm_step_err, norm_params = alt_mass_step_calc(mu=panplus_data[:, hdr.index('MU')].astype(float),
                                                               mu_err=panplus_data[:, hdr.index('MUERR')].astype(float),
                                                               mass=panplus_data[:, hdr.index('HOST_LOGMASS')].astype(
                                                                   float),
                                                               z=panplus_data[:, hdr.index('zCMB')].astype(float),
                                                               cut=10, quiet=False)

    print('\nMass Step of 91bg-like HiloCATS - snpy')
    snpy_step, snpy_step_err, snpy_params = alt_mass_step_calc(mu=hicat_snpy_data[:, 8].astype(float),
                                                               mu_err=hicat_snpy_data[:, 9].astype(float),
                                                               mass=hicat_snpy_data[:, 16].astype(float),
                                                               z=hicat_snpy_data[:, 4].astype(float),
                                                               cut=10, quiet=False)

    print('\nMass Step of 91bg-like HiloCATS - salt')
    salt_step, salt_step_err, salt_params = alt_mass_step_calc(mu=hicat_salt_data[:, 16].astype(float),
                                                               mu_err=hicat_salt_data[:, 17].astype(float),
                                                               mass=hicat_salt_data[:, 18].astype(float),
                                                               z=hicat_salt_data[:, 4].astype(float),
                                                               cut=10, quiet=False)

    print('\nMass Step of 91bg-like HiloCATS - merged')
    merged_step, merged_step_err, merged_params = alt_mass_step_calc(mu=hicat_merged_data[:, 3].astype(float),
                                                                     mu_err=hicat_merged_data[:, 4].astype(float),
                                                                     mass=hicat_merged_data[:, 5].astype(float),
                                                                     z=hicat_merged_data[:, 1].astype(float),
                                                                     cut=10, quiet=False)

    # Plot mass_steps
    plt.figure(figsize=(9, 5))

    plt.errorbar(x=snpy_params['scatter'], y=norm_step - snpy_step, yerr=snpy_step_err, fmt='o')
    plt.text(snpy_params['scatter'], norm_step - snpy_step, 'snpy', size='x-small', va='bottom')

    plt.errorbar(x=salt_params['scatter'], y=norm_step - salt_step, yerr=salt_step_err, fmt='o')
    plt.text(salt_params['scatter'], norm_step - salt_step, 'salt', size='x-small', va='bottom')

    plt.errorbar(x=merged_params['scatter'], y=norm_step - merged_step, yerr=merged_step_err, fmt='o')
    plt.text(merged_params['scatter'], norm_step - merged_step, 'merged', size='x-small', va='bottom')

    plt.title('Mass Steps of SNooPy & SALT2 v. Normal SNe Ia Mass Step from Pantheon+')
    plt.xlabel('Scatter of Hubble Residuals');
    plt.ylabel('91bg Mass Step - Norm Ia Mass Step')
    # plt.axhline(y=0)
    plt.ylim(-0.45, 0.45)
    # plt.ylim(-0.15, 0.15)
    plt.xlim(0, 0.9)
    plt.show()
    return
# ------------------------------------------------------------------------------------------------------------------- #
def convert_panth_to_new_type():
    panplus_data = np.genfromtxt('txts/hlsp_ps1cosmo_panstarrs_gpc1_all_model_v1_ancillary-g10.fitres.txt', skip_header=69, dtype=str)[:,1:]  # Removes 'SN:' at begining of each line
    # hdr = ('CID', 'CIDint', 'IDSURVEY', 'TYPE', 'FIELD', 'CUTFLAG_SNANA', 'zHEL', 'zHELERR', 'zCMB', 'zCMBERR', 'zHD',
    #        'zHDERR', 'VPEC', 'VPECERR', 'MWEBV', 'HOST_LOGMASS', 'HOST_LOGMASS_ERR', 'HOST_sSFR', 'HOST_sSFR_ERR',
    #        'PKMJDINI', 'SNRMAX1', 'SNRMAX2', 'SNRMAX3', 'PKMJD', 'PKMJDERR', 'x1', 'x1ERR', 'c', 'cERR', 'mB', 'mBERR',
    #        'x0', 'x0ERR', 'COV_x1_c', 'COV_x1_x0', 'COV_c_x0', 'NDOF', 'FITCHI2', 'FITPROB', 'RA', 'DEC', 'HOST_RA',
    #        'HOST_DEC', 'HOST_ANGSEP', 'TGAPMAX', 'TrestMIN', 'TrestMAX', 'ELU', 'HOSTGAL_SFR', 'HOSTGAL_SFR_ERR',
    #        'HOSTGAL_sSFR', 'HOSTGAL_sSFR_ERR', 'CUTMASK', 'MU', 'MUMODEL', 'MUERR', 'MUERR_RENORM', 'MUERR_RAW',
    #        'MUERR_VPEC', 'MURES', 'MUPULL', 'M0DIF', 'M0DIFERR', 'CHI2', 'biasCor_mu', 'biasCorErr_mu',
    #        'biasCor_mB', 'biasCor_x1', 'biasCor_c', 'biasScale_muCOV', 'IDSAMPLE', 'IZBIN')
    hdr = ('CID CIDint IDSURVEY TYPE FIELD CUTFLAG_SNANA zCMB zCMBERR zHD zHDERR VPEC VPEC_ERR HOST_LOGMASS '
           'HOST_LOGMASS_ERR SNRMAX1 SNRMAX2 SNRMAX3 PKMJD PKMJDERR x1 x1ERR c cERR mB mBERR x0 x0ERR COV_x1_c COV_x1_x0 '
           'COV_c_x0 NDOF FITCHI2 FITPROB RA DECL TGAPMAX TrestMIN TrestMAX MWEBV MU MUMODEL MUERR MUERR_RAW MURES '
           'MUPULL ERRCODE biasCor_mu biasCorErr_mu biasCor_mB biasCor_x1 biasCor_c biasScale_muCOV IDSAMPLE ').split(' ')
    panplus_data = panplus_data[panplus_data[:, hdr.index('zCMB')].astype(float) < 0.08]  # Only low-z SNe

    with open('txts/PanPlus_adjusted.txt', 'w') as f:
        new_hdr = ['zCMB', 'zCMBERR', 'HOST_LOGMASS', 'HOST_LOGMASS_ERR', 'x1', 'x1ERR', 'c', 'cERR', 'x0', 'x0ERR', 'RA', 'DECL', 'MU', 'MUERR']
        f.write('objname, origin, z_cmb, z_cmb_err, hostMass, hostMass_err, x1, x1_err, c, c_err, x0_err, x0_err, ra, dec, mu, mu_err\n')
        for i in range(len(panplus_data)):
            line = panplus_data[i, 0] + ', Pan+'
            for cat in new_hdr:
                line += ', '+str(panplus_data[i, hdr.index(cat)].astype(float))
            line += '\n'
            f.write(line)
    return
def convert_dr3_to_new_type():
    dr3_data = np.genfromtxt('txts/DR3_fits.dat', skip_header=1, dtype=str)
    hdr = ['SN', 'st', 'e_st', 'zhelio', 'zcmb', 'Tmax', 'eTmax', 'Bmax', 'e_Bmax', 'Vmax', 'e_Vmax', 'umax', 'e_umax',
           'gmax', 'e_gmax', 'rmax', 'e_rmax', 'imax', 'e_imax', 'Ymax', 'e_Ymax', 'Jmax', 'e_Jmax', 'Hmax', 'e_Hmax',
           'EBVhost', 'e_EBVhost', 'Rv', 'e_Rv', 'C(Rv)(EBVhost)', 'c0', 'c1', 'c2', 'c3', 'c4', 'c5']
    with open('txts/dr3_adjusted.txt', 'w') as f:
        new_hdr = ['zhelio', 'zcmb', 'st', 'e_st', 'EBVhost', 'e_EBVhost', 'Tmax', 'eTmax']
        f.write('objname, origin, z, z_cmb, st, st_err, EBVhost, EBVhost_err, Tmax, Tmax_err\n')
        for i in range(len(dr3_data)):
            line = dr3_data[i, 0] + ', dr3'
            for cat in new_hdr:
                line += ', '+str(dr3_data[i, hdr.index(cat)].astype(float))
            line += '\n'
            f.write(line)
    return
def covert_normIa_to_my_type():
    normIa_data = np.genfromtxt('txts/txts/normIaparams.txt', skip_header=1, dtype=str)
    return
# ------------------------------------------------------------------------------------------------------------------- #
def help():
    """
    Call to get examples of calls for smart_fit()
    """
    print('===========================================================================================================')
    print('Ex. Indivisual:', "smart_fit(fit_type='indv', data_set='CSP', algo='snpy', path='data/CSP/SN2008bi_snpy.txt')")
    print('------')
    print('Ex. Batch:', "smart_fit(fit_type='batch', 'CSP', algo='snpy'")
    print('------')
    print('Ex. Combined:', "smart_fit(fit_type='combined', algo='snpy')")
    print('===========================================================================================================')
    return
def smart_fit(fit_type, data_set, algo, path=None, save_loc='', dmag_max=0.00, dflux_max=0.00):
    """
    An function to easily fit data using both algorithms.
    :param fit_type: [str] type of fitting protocol to use
    :param data_set: [str] name of data set
    :param algo: [str] algorithms to fit data with (snpy/salt)
    :param path: [str] Default: None, data location for indivisual fits (only call for fit_type='indv')
    :param save_loc: [str] location to save paramater data of SNe fits
    :param dmag_max: [float] Default: 0, magnitude error cut off for taking in data
    :param dflux_max: [float] Default: 0,  flux error cut off for taking in data
    :return: SN [list], Array of sn91bg() classes from fitting call. Returns list of 'None' if unsuccessful.
    """
    fit_type = fit_type.lower()
    if fit_type == 'indv':
        tempSN = sn91bg().make_class(data_set=data_set, path=path, dmag_max=dmag_max, dflux_max=dflux_max)
        tempSN.fit(algo)
        SNe = [tempSN]
        if SNe[0] == None:
            print('[!!!] Fit failed!')
            return [None]
    elif fit_type == 'batch':
        SNe = batch_fit(data_set, algo=algo,
                        dmag_max=dmag_max, dflux_max=dflux_max)
    elif fit_type == 'combiend':
        SNe = combined_fit(algo=algo,
                           dmag_max=dmag_max, dflux_max=dflux_max)
    else:
        raise ValueError(
            "Invalid type selected ['indv'/'batch'/'combined']")

    # Saving
    if fit_type != 'indv' and len(SNe) != 0:
        if len(save_loc) == 0:
            save_loc = 'output/'+fit_type+'_'+data_set.lower()+'_'+algo+'_'+'params.txt'
        save_params_to_file(save_loc, SNe) # Save params
        sample_cutter(save_loc, algo) # Cut sample
    return SNe

if __name__ == '__main__':
    start = systime.time() # Runtime tracker

    smart_fit('indv', 'CSP', 'snpy', path='data/CSP/SN2005ke_snpy.txt')

    # mass_step_residual_plot('output/merged_params_cut.txt',
    #                         title="Hubble Residuals vs. Host Stellar Mass of CSP-ATLAS-ZTF 91bg-like SNe Ia (salt-snpy)\n")
    # mass_step_residual_plot('HiCAT_DR3_params.txt',
    #                         title="Hubble Residuals vs. Host Stellar Mass of CSP Normal SNe Ia\n")



    # algo = 'snpy'
    # SNe, files = [], glob.glob('data/CSPdata/*.txt')
    # for path in files:
    #     print('[', files.index(path) + 1, '/', len(files), ']')
    #     print('-----------------------------------------------------------------------------------------------')
    #     tempSN = sn91bg().make_class(data_set='CSP', path=path)
    #     tempSN.fit(algo)
    #
    #     # Check if fit failed
    #     if (tempSN != None and
    #             tempSN.params['mu']['value'] > 0 and
    #             'hostMass' in tempSN.params and
    #             tempSN.params['hostMass']['value'] > 0):
    #         SNe.append(tempSN)
    #     # break
    # print('Sucessfully fit [', len(SNe), '/', len(files), ']!')
    #
    # # Save params to file
    # save_params_to_file('HiCAT_DR3_params.txt', SNe)
    #
    # # Cut sample
    # sample_cutter('HiCAT_DR3_params.txt', algo=algo)

    # smart_fit('COMBINED', 'COMBINED', 'salt', dmag_max=1)
    #
    # merge_snpy_salt('output/COMBINED_combined_snpy_params_cut.txt',
    #                 'output/COMBINED_combined_salt_params_cut.txt',
    #                 'output/merged_params_cut.txt')
    #
    # make_display_plots()

    # hicat_data = np.genfromtxt('output/COMBINED_combined_salt_params_cut.txt', delimiter=', ', skip_header=1, dtype=str)
    # panplus_data = np.genfromtxt('txts/PanPlus_Latest.FITRES', skip_header=69, dtype=str)[:, 1:] # Removes 'SN:' at begining of each line
    # hdr = ('CID CIDint IDSURVEY TYPE FIELD CUTFLAG_SNANA zHEL zHELERR zCMB zCMBERR zHD zHDERR VPEC VPECERR MWEBV '
    #        'HOST_LOGMASS HOST_LOGMASS_ERR HOST_sSFR HOST_sSFR_ERR PKMJDINI SNRMAX1 SNRMAX2 SNRMAX3 PKMJD PKMJDERR '
    #        'x1 x1ERR c cERR mB mBERR x0 x0ERR COV_x1_c COV_x1_x0 COV_c_x0 NDOF FITCHI2 FITPROB RA DEC HOST_RA HOST_DEC '
    #        'HOST_ANGSEP TGAPMAX TrestMIN TrestMAX ELU HOSTGAL_SFR HOSTGAL_SFR_ERR HOSTGAL_sSFR HOSTGAL_sSFR_ERR CUTMASK '
    #        'MU MUMODEL MUERR MUERR_RENORM MUERR_RAW MUERR_VPEC MURES MUPULL M0DIF M0DIFERR CHI2 biasCor_mu biasCorErr_mu '
    #        'biasCor_mB biasCor_x1 biasCor_c biasScale_muCOV IDSAMPLE IZBIN').split(' ')
    # panplus_data = panplus_data[panplus_data[:, hdr.index('zCMB')].astype(float) < np.max(hicat_data[:, 4].astype(float))] # Only low-z SNe
    #
    # param_hist_compare(panplus_data[:, hdr.index('x0')].astype(float), hicat_data[:, 10].astype(float), [0.0005, 0.001],
    #                    xrange=[0, 0.02], title='Normal-SNe-Ia v. 91bg-like -- x0')
    # param_hist_compare(panplus_data[:, hdr.index('x1')].astype(float), hicat_data[:, 12].astype(float), [0.125, 0.4],
    #                    xrange=None, title='Normal-SNe-Ia v. 91bg-like -- x1')
    # param_hist_compare(panplus_data[:, hdr.index('c')].astype(float), hicat_data[:, 14].astype(float), [0.008, 0.2],
    #                    xrange=None, title='Normal-SNe-Ia v. 91bg-like -- c')
    #


    # merge_snpy_salt('output/COMBINED_combined_snpy_params_cut.txt',
    #                 'output/COMBINED_combined_salt_params_cut.txt',
    #                 'mass_step_scatter.txt')

    # mass_step_calc('mass_step_scatter.txt', cut=10, quiet=False)
    # mass_step_calc('mass_step_scatter.txt', cut='median', quiet=False)

    # hicat_data = np.genfromtxt('output/COMBINED_combined_snpy_params_cut.txt', delimiter=', ', skip_header=1, dtype=str)
    # dr3_data = np.genfromtxt('txts/DR3_fits.dat', dtype=str, skip_header=1)
    # param_hist_compare(hicat_data[:, 14].astype(float), dr3_data[:, 25].astype(float), 0.035,
    #                    xrange=None, title='91bg-like v. Normal-SNe-Ia -- EBVhost')

    # sample_cutter('output/COMBINED_combined_salt_params.txt', algo='salt')
    # residual_plotter('output/COMBINED_combined_salt_params_cut.txt', x_params='z', algo='snpy', labels=True)

    # data = np.genfromtxt('output/COMBINED_combined_snpy_params.txt', delimiter=', ', skip_header=1, dtype=str)
    # low_mass_data = data[data[:, 16].astype(float) < 10]
    # for i in range(len(low_mass_data)):
    #     print(low_mass_data[i, 0], '|', low_mass_data[i, 10], '|', low_mass_data[i, 14])
    # print('-------')
    # print('Stetch:', np.average(low_mass_data[:, 10].astype(float)),
    #       np.min(low_mass_data[:, 10].astype(float)), np.max(low_mass_data[:, 10].astype(float)))
    # print('EBVhost:', np.average(low_mass_data[:, 14].astype(float)),
    #       np.min(low_mass_data[:, 14].astype(float)), np.max(low_mass_data[:, 14].astype(float)))


    # SNe = []
    # for n in ['2007ba']:
    #     SNe.append(sn91bg().load_from_file('saved/snpy/combined/classes/'+n+'_class.txt'))
    # for SN in SNe:
    #     print(SN.params['EBVhost'], SN.params['st'], SN.params['Tmax'])

    # smart_fit('COMBINED', 'COMBINED', 'snpy', dmag_max=1)
    # smart_fit('COMBINED', 'COMBINED', 'salt', dmag_max=1)

    # sample_cutter('saved/snpy/combined/combined_params.txt', algo='snpy')
    # residual_plotter('saved/snpy/combined/combined_params_cut.txt', x_params='mass', algo='snpy', labels=True)

    # SN = sn91bg().load_from_file('saved/snpy/combined/classes/2021twa_class.txt')
    # SN.plot()

    # for n_cat in ['t0_err', 'x0', 'x0_err', 'x1', 'x1_err', 'c', 'c_err']:
    #     path = 'output/COMBINED_combined_salt_params.txt'
    #     data = np.genfromtxt(path, delimiter=', ', skip_header=1, dtype=str)
    #     with open(path, 'r') as f:
    #         hdr = f.readline().split(', ')
    #         hdr[-1] = hdr[-1][:-1]
    #
    #     cat = data[:, hdr.index(n_cat)].astype(float)
    #     ebv_avg, ebv_std = np.average(cat), np.std(cat)
    #     z_scores = (cat - ebv_avg) / ebv_std
    #
    #     print(n_cat, np.min(cat[z_scores < 3]), '< x <', np.max(cat[z_scores < 3]))
    #     # break

    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')

