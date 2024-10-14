# M.D. Woods - 10/11/24
import warnings
# warnings.simplefilter("ignore", UserWarning)
import os
import sys

import numpy
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
from astroquery.sdss import SDSS
from astropy.stats import sigma_clip
from astropy.stats import sigma_clipped_stats

from scripts import general as gen
from scripts import get_vpec

CONSTANTS = gen.get_constants()

class SN91bg:
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
        return (self.objname + ' | ' + self.originalname + ' | ' + self.origin + ' | ' + str(self.coords) +
                ' | (z | z_cmb): (' + str(self.z) + ' | ' + str(round(self.z_cmb, 2)) + ')')
    def print_info(self):
        # Header
        print('---------------------------------------------------------------------------------------------')
        print('|+| ' + self.objname + ' |+|\n' + self.originalname + ' | ' + self.origin + ' | ' + str(self.coords) +
              ' | (z | z_cmb): (' + str(self.z) + ' | ' + str(round(self.z_cmb, 4)) + ')')
        print('---------------------------------------------------------------------------------------------')

        # Print array ranges
        print('\tFilters: \t' + str(np.unique(self.filters)))
        print('\tFlux: \t\t(' + str(round(np.min(self.flux), 3)) + ' -- ' + str(round(np.max(self.flux), 3)) + ')' +
              ' | ' + 'dFlux: (' + str(round(np.min(self.dflux), 3)) + ' -- ' + str(round(np.max(self.dflux), 3)) + ')')
        print('\tMag: \t\t(' + str(round(np.min(self.mag), 3)) + ' -- ' + str(round(np.max(self.mag), 3)) + ')' +
              ' | ' + 'dMag: (' + str(round(np.min(self.dmag), 3)) + ' -- ' + str(round(np.max(self.dmag), 3)) + ')')
        print('\tMJDs-MJDe: \t(' + str(self.period[0]) + ' -- ' + str(self.period[1]) + ')')
        print('---------------------------------------------------------------------------------------------')

        # Print params
        if len(self.params) > 0:
            for p in self.params:
                print('\t'+p+': ' + str(round(self.params[p]['value'], 3)) +
                      ' +/- ' + str(round(self.params[p]['err'], 3)))
        print('---------------------------------------------------------------------------------------------')
    def plot(self, y_type='mag', save_loc='', zoom=0, subplots=False, date_lines=True):
        print('[+++] Plotting LC of '+self.objname+'...')
        filter_dict = {'u': 'teal', 'g': 'green', 'r': 'red', 'i': 'indigo', 'B': 'blue',
                       'V0': 'violet', 'V1': 'purple', 'V': 'red', 'Y': 'goldenrod', 'Hdw': 'tomato', 'H': 'salmon',
                       'J': 'aquamarine', 'Jrc2': 'cadetblue', 'Jdw': 'turquoise', 'Ydw': 'olive',
                       'c': 'cyan', 'o': 'orange', 'ZTF_g': 'green', 'ZTF_r': 'red', 'ZTF_i': 'indigo'}

        # Plot
        unique_filters, num_plts = np.unique(self.filters), len(np.unique(self.filters))
        if not subplots:
            plt.figure(figsize=(10, 4), constrained_layout=True)
            size = None
        elif self.origin == 'CSP':
            plt.figure(figsize=(25, 8), constrained_layout=True)
            size = (2, 6)
        elif self.origin == 'ATLAS':
            plt.figure(figsize=(12, 4), constrained_layout=True)
            size = (1, 2)
        elif self.origin == 'ZTF':
            plt.figure(figsize=(25, 5), constrained_layout=True)
            size = (1, 3)
        elif self.origin == 'ATLAS-ZTF':
            plt.figure(figsize=(25, 8), constrained_layout=True)
            size = (2, 5)
        else:
            raise ValueError('[!!!] Origin not valid!'
                             )

        # Plot for each filter
        for i in range(num_plts):
            if size is not None:
                plt.subplot(size[0], size[1], i+1)
            indexes = np.where(self.filters == unique_filters[i])[0]
            if y_type == 'mag':
                plt.errorbar(self.time[indexes], self.mag[indexes], yerr=self.dmag[indexes], fmt='o', ms=4, elinewidth=0.3,
                             color=filter_dict[self.filters[indexes][0]], label=self.filters[indexes][0])
            elif y_type == 'flux':
                plt.errorbar(self.time[indexes], self.flux[indexes], yerr=self.dflux[indexes],
                             fmt='o', ms=4, elinewidth=0.3,
                             color=filter_dict[self.filters[indexes][0]], label=self.filters[indexes][0])

            # Format
            if size is not None and i > 0:
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
            plt.axvline(x=float(self.discovery_date),
                        color='peru', ls='--', label='Discovery Date', linewidth=3, alpha=0.5)

        if zoom > 0:
            if 'Tmax' in list(self.params.keys()):
                plt.xlim(self.params['Tmax']['value']-zoom, self.params['Tmax']['value']+zoom)
            elif 't0' in list(self.params.keys()):
                plt.xlim(self.params['t0']['value']-zoom, self.params['t0']['value']+zoom)

            plt.xlabel('MJD')
            plt.legend()
        if y_type == 'mag':
            plt.gca().invert_yaxis()
        elif y_type == 'flux':
            plt.ylim(0)
        plt.suptitle('Lightcurve -- '+self.objname+' | '+self.originalname+' -- '+y_type)
        if len(save_loc) != 0:
            print('[+++] Saving to '+save_loc)
            plt.savefig(save_loc)
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
                hdr[0], hdr[1], hdr[2],
                (float(hdr[3]), float(hdr[4])), float(hdr[5]), float(hdr[6]), float(hdr[7][:-1]))

            f.readline()

            line = f.readline()
            while '+++' not in line:
                line = line.split(', ')
                self.params.update({line[0]: {'value': float(line[1]), 'err': float(line[2])}})
                line = f.readline()

            self.zp = np.array(f.readline().split(', ')[1:]).astype(float)
            self.filters = np.array(f.readline().split(', ')[1:]).astype(str)
            self.filters[-1] = self.filters[-1][:-1]  # Removes the /n from the end of the line
            self.time = np.array(f.readline().split(', ')[1:]).astype(float)
            self.flux = np.array(f.readline().split(', ')[1:]).astype(float)
            self.dflux = np.array(f.readline().split(', ')[1:]).astype(float)
            self.mag = np.array(f.readline().split(', ')[1:]).astype(float)
            self.dmag = np.array(f.readline().split(', ')[1:]).astype(float)

            self.period = (np.min(self.time), np.max(self.time))
        return self
    # -----------------------------------------------------------------------------------------------------------------
    def fit(self, algo='snpy', save_loc=None):
        if save_loc is None:
            save_loc = CONSTANTS[algo+'_'+self.origin.lower() + '_saved_loc']

        # Choose algorithm
        if algo == "snpy":
            self.write_snpy_ascii(save_loc=save_loc + 'ascii/')
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
    def make_class(self, data_set, path, dmag_max=0.00, dflux_max=0.00):
        if data_set == 'CSP':
            print('[+++] Creating class using CSP data...')
            print(path)
            # Header elements
            with open(path, 'r') as f:
                objname, z, ra, dec = f.readline().split(' ')
                objname, z, ra, dec = objname[2:], float(z), float(ra), float(dec[:-1])
            originalname = path.split('/')[-1].split('_')[0]

            # Query TNS for transient details
            objname, z_void, discdate = gen.TNS_details(ra, dec)  # Voiding the redshift from TNS, its None for some reason

            # Make class
            temp_sn = SN91bg(objname=objname,
                             originalname=originalname,
                             coords=(ra, dec),
                             z=z,
                             origin='CSP',
                             discovery_data=discdate)

            # Set arrays
            class_filter = ''
            with open(path, 'r') as f:
                f.readline()  # Skips header
                for line in f.readlines():
                    data_line = line[:-1].split(' ')
                    if len(data_line) == 2:
                        class_filter = data_line[1]
                    elif len(data_line) >= 3:
                        if len(data_line) == 4:
                            data_line = data_line[1:]
                        n_time, n_mag, n_dmag = float(data_line[0]), float(data_line[1]), float(data_line[2])
                        zp = float(CONSTANTS['csp_zpts_' + class_filter])
                        n_time = n_time + 53000  # JD to MJD

                        n_flux = 10 ** ((n_mag - zp) / -2.5)
                        n_dflux = np.abs(n_flux) * np.log(10) * ((1 / 2.5) * n_dmag)

                        temp_sn.zp = np.append(temp_sn.zp, zp)
                        temp_sn.filters = np.append(temp_sn.filters, class_filter)
                        temp_sn.time = np.append(temp_sn.time, n_time)
                        temp_sn.flux = np.append(temp_sn.flux, n_flux)
                        temp_sn.dflux = np.append(temp_sn.dflux, n_dflux)
                        temp_sn.mag = np.append(temp_sn.mag, n_mag)
                        temp_sn.dmag = np.append(temp_sn.dmag, n_dmag)
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
            temp_sn = SN91bg(objname=objname,
                             originalname=path.split('/')[-1].split('.')[0],
                             coords=(ra, dec),
                             z=np.nan if z == 'None' else float(z),
                             origin='ATLAS',
                             discovery_data=discdate)

            # Set arrays
            temp_sn.zp = data[:, 7].astype(float)
            temp_sn.filters = data[:, 6]
            temp_sn.time = data[:, 8]
            temp_sn.flux = data[:, 16]
            temp_sn.dflux = data[:, 17]
            temp_sn.mag = data[:, 3]
            temp_sn.dmag = data[:, 4]
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
            temp_sn = SN91bg(objname=objname,
                             originalname=originalname,
                             coords=(ra, dec),
                             z=z,
                             origin='ZTF',
                             discovery_data=discdate)

            # Set arrays
            temp_sn.zp = zp
            temp_sn.filters = filters
            temp_sn.time = time
            temp_sn.flux = flux
            temp_sn.dflux = dflux
            temp_sn.mag = mag
            temp_sn.dmag = dmag
        else:
            raise ValueError("Data set '" +
                             data_set + "' not recognized")

        # Clean data
        temp_sn.clean_data(dmag_max, dflux_max)

        # Final stage
        if temp_sn.period is None:
            print('[!!!] No valid points found in file!')
            return None
        else:
            print('      Class created successfully!')
            return temp_sn
    def clean_data(self, dmag_max=0.00, dflux_max=0.00):
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
            n_zp, n_filters, n_time = self.zp[n], self.filters[n], self.time[n]
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
    def write_snpy_ascii(self, save_loc='default/'):
        print('[+++] '+self.objname+' -- Saving data to ASCII files for SNooPy...')
        filter_dict = {'o': 'ATri', 'c': 'ATgr',
                       'ZTF_g': 'g', 'ZTF_r': 'r', 'ZTF_i': 'i',
                       'B': 'B', 'H': 'H', 'J': 'J', 'Jrc2': 'Jrc2', 'V': 'V', 'V0': 'V0', 'Y': 'Y', 'Ydw': 'Ydw', 'g': 'g', 'i': 'i', 'r': 'r', 'u': 'u'}
        with open(save_loc + self.objname + '_snpy.txt', 'w') as f:
            # Line 1 -- Objname, Helio-Z, RA, Dec (Ex. SN1981D 0.005871 50.65992 -37.23272)
            f.write(str(self.objname)+' '+str(self.z)+' '+str(self.coords[0])+' '+str(self.coords[1])+'\n')
            for f_w in np.unique(self.filters):
                f_indexes = np.where(self.filters == f_w)[0]
                f.write('filter ' + filter_dict[f_w] + '\n')
                for i in f_indexes:
                    # filter photometry block -- Date (JD/MJD), mag, err (i.e. 674.8593 12.94 0.11)
                    f.write(str(self.time[i]) + '\t' + str(self.mag[i]) + '\t' + str(self.dmag[i]) + '\n')
        print('      Saved file to '+save_loc + self.objname + '_snpy.txt')
        return
    def snpy_fit(self, save_loc, use_saved=False, show_plot=True, quiet=False):
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
        except Exception as error:
            self.params.update({'mu': {'value': 0.00, 'err': 0.00}})
            print('[!!!] Failed to load ASCII file -- ', error)
            return
        # n_s.k_version = '91bg'
        n_s.choose_model('EBV_model2', stype='st')
        n_s.set_restbands()  # Auto pick appropriate rest-bands

        # Remove empty filters -- fix for 'ValueError: attempt to get argmin of an empty sequence'
        for class_filter in list(n_s.data.keys()):
            if len(n_s.data[class_filter].magnitude) == 0:
                del n_s.data[class_filter]
            elif self.origin == 'CSP' and class_filter in ['u', 'Y', 'J', 'H', 'Jrc2', 'Ydw']:
                print('[***] Special Process for CSP! Removing ' + class_filter + '...')
                del n_s.data[class_filter]
        print('      Best filters:', list(n_s.data.keys()))

        for i in range(5):
            try:
                if self.origin == 'CSP':
                    initial_filters = []
                    for fil in ['B', 'V', 'g']:
                        if fil in list(n_s.data.keys()):
                            initial_filters.append(fil)
                    print('[***] Special Process for CSP! Fitting as '+str(initial_filters)+' -> remaining...')

                    n_s.fit(initial_filters, dokcorr=True, k_stretch=False, reset_kcorrs=True,
                            **{'mangle': 1, 'calibration': 0})
                    n_s.fit(bands=None, dokcorr=True, k_stretch=False, reset_kcorrs=True,
                            **{'mangle': 1, 'calibration': 0})
                else:
                    n_s.fit(bands=None, dokcorr=True, k_stretch=False, reset_kcorrs=True,
                            **{'mangle': 1, 'calibration': 0})
                n_s.save(save_path)

                # Save parameters
                for j in range(len(param_names)):
                    self.params.update({param_names[j]: {'value': n_s.parameters[snpy_param_names[j]],
                                                         'err': n_s.errors[snpy_param_names[j]]}})

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
            # details = gen.TNS_details(self.coords[0], self.coords[1])
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
        local_coords = SkyCoord(self.coords[0], self.coords[1], unit="deg")
        galac_coords = local_coords.transform_to(Galactic())

        # Get CMB redshift
        helio_corr = (float(CONSTANTS['cmb_v_helio']) / float(CONSTANTS['cmb_c']) *
                      ((np.sin(galac_coords.b.deg) * np.sin(float(CONSTANTS['cmb_b_h'])) + np.cos(galac_coords.b.deg) *
                        np.cos(float(CONSTANTS['cmb_b_h'])) * np.cos(galac_coords.l.deg - float(CONSTANTS['cmb_l_h'])))))
        corr_term = 1 - helio_corr
        self.z_cmb = (1 + self.z) / corr_term - 1

        # Peculiar Velocity Correction -- using 'get_vpec.py' from David
        VP = get_vpec.VelocityCorrection(f"twomass++_velocity_LH11.npy")
        self.z_cmb = VP.correct_redshift(self.z, 0, local_coords.galactic.l.deg, local_coords.galactic.b.deg)
        vpec, vpec_sys = get_vpec.main(self.coords[0], self.coords[1], self.z_cmb)
        self.z_cmb += vpec / 3e5

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
        gMag, iMag, iAbsMag = -999.00, -999.00, -999.00
        gMagErr, iMagErr, iAbsMagErr = -999.00, -999.00, -999.00
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
            if result is None:
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

        # Mass Calculation -- Taylor et al. 2011 -- eq. 8
        host_mass = (1.15 + (0.7 * (gMag - iMag)) - (0.4 * iAbsMag))

        # Error Propagation
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

# Fitting Functions ------------------------------------------------------------------------------------------------- #
def norm_fit(algo='snpy', save_loc=None, dmag_max=0.00, dflux_max=0.00):
    SNe, files = [], glob.glob('data/CSPdata/*.txt')
    for path in files:
        print('[', files.index(path) + 1, '/', len(files), ']')
        print('-----------------------------------------------------------------------------------------------')
        tempSN = SN91bg().make_class(data_set='CSP', path=path, dmag_max=dmag_max, dflux_max=dflux_max)
        if tempSN is not None:
            tempSN.fit(algo)

        # Check if fit failed
        if (tempSN is None or
                tempSN.params['mu']['value'] <= 0 or
                'hostMass' not in tempSN.params or
                tempSN.params['hostMass']['value'] <= 0):
            print('[-----] Failed!')
        else:
            print('[+++++] Success!')
            SNe.append(tempSN)

    print('Sucessfully fit [', len(SNe), '/', len(files), ']!')
    if save_loc is not None:
        save_params_to_file(save_loc, SNe)
    return SNe
def combined_fit(algo='snpy', dmag_max=0.00, dflux_max=0.00):
    sys.stdout = open(os.devnull, 'w')  # Lots of unnecessary output

    # Load and set up files
    csp_files, atlas_files, ztf_files = glob.glob('data/CSP/*.txt'), glob.glob('data/ATLAS/*.txt'), glob.glob('data/ZTF/*.txt')
    CSP_SNe, ATLAS_SNe, ZTF_SNe, atlas_names = {}, {}, {}, []
    for file in csp_files:
        tempSN = SN91bg().make_class('CSP', file, dmag_max, dflux_max)
        if tempSN is not None:
            CSP_SNe.update({tempSN.objname: tempSN})
    for file in atlas_files:
        tempSN = SN91bg().make_class('ATLAS', file, dmag_max, dflux_max)
        if tempSN is not None:
            ATLAS_SNe.update({tempSN.objname: tempSN})
    for file in ztf_files:
        tempSN = SN91bg().make_class('ZTF', file, dmag_max, dflux_max)
        if tempSN is not None:
            ZTF_SNe.update({tempSN.objname: tempSN})

    sys.stdout = sys.__stdout__

    # List of overlap
    atlas_ztf_list = []
    for element in list(ATLAS_SNe.keys()):
        if element in list(ZTF_SNe.keys()):
            atlas_ztf_list.append(element)

    # List of all unique SNe
    combined_list = np.unique(np.hstack((np.hstack((list(CSP_SNe.keys()), list(ATLAS_SNe.keys()))),
                                         list(ZTF_SNe.keys()))))

    # Combine filters and data
    combined_SNe = []
    for name in combined_list:
        if name in atlas_ztf_list:
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
    for n_SN in combined_SNe:
        print('-----------------------------------------------------------------------------------------------------')
        print('[', combined_SNe.index(n_SN) + 1, '/', len(combined_SNe), '] Fitting data for ' + n_SN.objname + ' [' + algo + '] [' + n_SN.origin + ']...')
        n_SN.fit(algo=algo, save_loc=CONSTANTS[algo + '_combined_saved_loc'])

        # Check if fit failed
        if (n_SN is None or
                n_SN.params['mu']['value'] <= 0 or
                'hostMass' not in n_SN.params or
                n_SN.params['hostMass']['value'] <= 0):
            print('[-----] Failed!')
        else:
            print('[+++++] Success!')
            fit_combined_SNe.append(n_SN)

    print('=====================================================================================================\n',
          'Successfully fit [', len(fit_combined_SNe), '/', len(combined_SNe), ']!\n',
          '=====================================================================================================\n')
    return fit_combined_SNe
def batch_fit(data_set, algo='snpy', dmag_max=0.00, dflux_max=0.00):
    SNe, files = [], glob.glob(CONSTANTS[data_set.lower()+'_data_loc'] + '*.txt')
    for path in files:
        print('[', files.index(path) + 1, '/', len(files), ']')
        print('-----------------------------------------------------------------------------------------------')
        tempSN = SN91bg().make_class(data_set=data_set, path=path, dmag_max=dmag_max, dflux_max=dflux_max)
        if tempSN is not None:
            tempSN.fit(algo)

        # Check if fit failed
        if (tempSN is None or
                tempSN.params['mu']['value'] <= 0 or
                'hostMass' not in tempSN.params or
                tempSN.params['hostMass']['value'] <= 0):
            print('[-----] Failed!')
        else:
            print('[+++++] Success!')
            SNe.append(tempSN)

    print('Sucessfully fit [', len(SNe), '/', len(files), ']!')
    return SNe
def batch_load(data_set, algo='snpy'):
    SNe = []
    for path in glob.glob('saved/' + algo + '/' + data_set.lower() + '/classes/*_class.txt'):
        SNe.append(SN91bg().load_from_file(path))
    return SNe

# File Saving / Augmentation ---------------------------------------------------------------------------------------- #
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
        for n_SN in SNe:
            line = (n_SN.objname + ', ' + str(n_SN.coords[0]) + ', ' + str(n_SN.coords[1]) + ', ' +
                    str(n_SN.z) + ', ' + str(n_SN.z_cmb) + ', ' +
                    str(n_SN.period[0]) + ', ' + str(n_SN.period[1]) + ', ' + str(n_SN.origin))
            for param in n_SN.params:
                line += ', ' + str(n_SN.params[param]['value']) + ', ' + str(n_SN.params[param]['err'])
            f.write(line + '\n')
    return
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
        cuts = {'z': 0.015, 'EBVhost': (-0.2, 0.3), 'EBVhost_err': 0.1, 'st': (-999, 0.8), 'st_err': 0.1, 'Tmax_err': 1, 'mu_err': 0.2}
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
                    (data[:, hdr.index('Tmax_err')].astype(float) < cuts['Tmax_err']) &
                    (data[:, hdr.index('mu_err')].astype(float) < cuts['mu_err'])]
    elif algo == 'salt':
        print('[+++] Cutting sample for SALT data...')
        # cuts = {'z': 0.015, 'c': (-0.3, 1.0), 'c_err': 0.1, 'x1_err': 1.0, 't0_err': 1, 'mu_err': 0.2}
        cuts = {'z': 0.015, 'c': (-0.3, 0.6), 'c_err': 0.1, 'x1_err': 0.5, 't0_err': 1, 'mu_err': 0.2}
        f_out = '      | '
        for c in cuts:
            f_out += c + ': ' + str(cuts[c]) + ' | '
        print(f_out)
        data = data[(data[:, hdr.index('z_cmb')].astype(float) > cuts['z']) &
                    (data[:, hdr.index('c')].astype(float) > cuts['c'][0]) &
                    (data[:, hdr.index('c')].astype(float) < cuts['c'][1]) &
                    (data[:, hdr.index('c_err')].astype(float) < cuts['c_err']) &
                    (data[:, hdr.index('x1_err')].astype(float) < cuts['x1_err']) &
                    (data[:, hdr.index('t0_err')].astype(float) < cuts['t0_err']) &
                    (data[:, hdr.index('mu_err')].astype(float) < cuts['mu_err'])]

    # Save to file
    with open(path[:-4]+'_cut.txt', 'w') as f:
        f_out = hdr[0]
        for h in hdr[1:]:
            f_out += ', ' + h
        f.write(f_out + '\n')
        for line in data[:]:
            f_out = line[0]
            for item in line[1:]:
                f_out += ', ' + str(item)
            f.write(f_out + '\n')
    print('      Cut file saved to...', path[:-4]+'_cut.txt')
    print('      [ '+str(len(data[:, 0]))+' / '+original_num+' ]')

    # Display Residual Scatter
    resid_scatter = sigma_clipped_stats(data[:, hdr.index('mu')].astype(float) -
                                   gen.current_cosmo().distmod(data[:, hdr.index('z_cmb')].astype(float)).value)[2]
    print('Hubble Residual Scatter:', resid_scatter)

    return
def merge_snpy_salt_params(snpy_path, salt_path, save_loc):
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

# Plotting Functions ------------------------------------------------------------------------------------------------ #
def resid_v_z(path, title='', save_loc=None):
    color_wheel = {'ZTF': '#D8973C', 'ATLAS': '#BD632F', 'CSP': '#72513B', 'ATLAS-ZTF': '#273E47',
                   'ZTF_SNPY': '#D8973C', 'ATLAS_SNPY': '#BD632F', 'CSP_SNPY': '#72513B', 'ATLAS-ZTF_SNPY': '#273E47',
                   'ZTF_SALT': '#D8973C', 'ATLAS_SALT': '#BD632F', 'CSP_SALT': '#72513B', 'ATLAS-ZTF_SALT': '#273E47',
                   'Pan+': 'BD632F', 'Histogram': '#3B5058',
                   '10dexfill': '#A4243B', 'mediandexfill': '#D8C99B'}
    fig, axs = plt.subplots(1, 2, figsize=(15, 6), gridspec_kw={'width_ratios': [10, 1]}, constrained_layout=True)

    # Pull data from saved text & header
    data = np.genfromtxt(path, delimiter=', ', skip_header=1, dtype=str)
    with open(path, 'r') as f:
        hdr = f.readline().split(', ')
        hdr[-1] = hdr[-1][:-1]

    # Set Arrays
    z = data[:, hdr.index('z_cmb')].astype(float)
    mass, mass_err = data[:, hdr.index('hostMass')].astype(float), data[:, hdr.index('hostMass_err')].astype(float)
    mu, mu_err = data[:, hdr.index('mu')].astype(float), data[:, hdr.index('mu_err')].astype(float)
    resid_mu, resid_mu_err = sigma_clip(mu - gen.current_cosmo().distmod(z).value, sigma=3.0), np.copy(mu_err)

    # Make main plot
    for origin in np.unique(data[:, hdr.index('origin')]):
        format_dict = {'marker': 'o', 'fmt': 'o', 'label': origin, 'alpha': 1, 'ms': 6}
        if 'SNPY' in origin:
            format_dict['label'] = origin[:-5]
        elif 'SALT' in origin:
            format_dict['label'], format_dict['marker'] = None, '^'

        indexes = np.where(data[:, hdr.index('origin')] == origin)[0]
        axs[0].errorbar(z[indexes], resid_mu[indexes], yerr=resid_mu_err[indexes],
                        color=color_wheel[origin], elinewidth=0.8, **format_dict)

    # Make histogram
    axs[1].hist(resid_mu, bins=int((np.max(resid_mu) - np.min(resid_mu)) / 0.1),  # Bin Width = 0.1
                orientation="horizontal", color=color_wheel['Histogram'])

    # Extra Info
    extra_info = '$\sigma$: '+str(round(np.std(resid_mu), 4)) + ', $n$: ' + str(len(resid_mu))
    if 'merged' in path:
        extra_info += r' | SALT3: $\triangle$, SNooPy: $\bigcirc$'
    axs[0].text(np.min(z), np.max(resid_mu),
                extra_info,
                horizontalalignment='left', verticalalignment='bottom')

    # Formatting
    fig.suptitle(title)
    axs[0].set(xlabel='Host Galaxy CMB Redshift', ylabel='Hubble Residuals (mag)')  # Sub-plot Labels
    axs[1].get_yaxis().set_visible(False)  # Turn off y-axis labels
    axs[0].legend(loc='best')

    # Saving Figure
    if save_loc is not None:
        print('Saved figure to... ', save_loc)
        plt.savefig(save_loc)
    plt.show()
    return
def resid_v_mass(path, cut=10, title='', labels=False, save_loc=None):
    # Pull data from saved text
    data = np.genfromtxt(path, delimiter=', ', skip_header=1, dtype=str)

    # Get header
    with open(path, 'r') as f:
        hdr = f.readline().split(', ')
        hdr[-1] = hdr[-1][:-1]

    fig, axs = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'width_ratios': [10, 1]}, constrained_layout=True)
    color_wheel = {'ZTF': '#81ADC8', 'ATLAS': '#EEDFAA', 'CSP': '#CD4631', 'ATLAS-ZTF': '#DEA47E',
                   'ZTF_SNPY': '#dcb160', 'ATLAS_SNPY': '#f4ac45', 'CSP_SNPY': '#af7b3f', 'ATLAS-ZTF_SNPY': '#694a38',
                   'ZTF_SALT': '#dcb160', 'ATLAS_SALT': '#f4ac45', 'CSP_SALT': '#af7b3f', 'ATLAS-ZTF_SALT': '#694a38',
                   'Pan+': 'maroon', 'Histogram': '#775A4A'}

    # Set Arrays
    mass, mass_err = data[:, hdr.index('hostMass')].astype(float), data[:, hdr.index('hostMass_err')].astype(float)
    mu, mu_err = data[:, hdr.index('mu')].astype(float), data[:, hdr.index('mu_err')].astype(float)
    z = data[:, hdr.index('z_cmb')].astype(float)
    resid_mu = sigma_clip(mu - gen.current_cosmo().distmod(z).value, sigma=3.0)
    resid_mu_err = np.copy(mu_err)

    # Make main plot
    for origin in np.unique(data[:, hdr.index('origin')]):
        if 'SALT' in origin:
            mk = 'x'
        else:
            mk = 'o'
        indexes = np.where(data[:, hdr.index('origin')] == origin)[0]

        axs[0].errorbar(mass[indexes], resid_mu[indexes], xerr=mass_err[indexes], yerr=resid_mu_err[indexes],
                        color=color_wheel[origin], fmt='o', label=origin, alpha=1, marker=mk)

        # Labels
        if labels:
            for i in range(len(resid_mu[indexes])):
                axs[0].text(mass[indexes][i], resid_mu[indexes][i], str(data[:, 0][indexes][i]), size='x-small', va='top')

    # Make histogram
    axs[1].hist(resid_mu, bins=30, orientation="horizontal", color=color_wheel['Histogram'])

    # Get Mass Step
    if cut == 'median':
        cut = round(np.median(mass), 4)
    mass_step_dict, resid_dict = mass_step_calc(mu, mu_err, mass, z, cut=cut)

    # Plot Mass Step
    plt_details = {'linestyle': '--', 'linewidth': 1.5, 'color': 'k'}
    fill_details = {'color': 'k', 'alpha': 0.2}
    axs[0].hlines(y=resid_dict['lower_resid']['value'], xmin=np.min(mass)-0.3, xmax=cut, **plt_details)  # Left
    axs[0].hlines(y=resid_dict['upper_resid']['value'], xmin=cut, xmax=np.max(mass)+0.3, **plt_details)  # Right
    axs[0].fill_between([np.min(mass)-0.3, cut], resid_dict['lower_resid']['value'] - resid_dict['lower_resid']['err'],
                        resid_dict['lower_resid']['value'] + resid_dict['lower_resid']['err'], **fill_details)  # Left
    axs[0].fill_between([cut, np.max(mass)+0.3], resid_dict['upper_resid']['value'] - resid_dict['upper_resid']['err'],
                        resid_dict['upper_resid']['value'] + resid_dict['upper_resid']['err'], **fill_details)  # Right
    axs[0].vlines(x=cut, ymin=resid_dict['lower_resid']['value'], ymax=resid_dict['upper_resid']['value'], **plt_details)

    # Formatting
    fig.suptitle(title +
                 'Scatter: ' + str(round(sigma_clipped_stats(resid_mu)[2], 2)) + ' | # of pts: ' + str(len(resid_mu)) +
                 ' | SNR: ' + str(round(np.sqrt(np.abs(np.average(resid_mu)) / np.abs(np.std(resid_mu_err))), 2)) +
                 '\nMass Step ('+str(cut)+'): ' + str(round(mass_step_dict['value'], 4)) + ' +/- ' + str(round(mass_step_dict['err'], 6)),
                 size='medium')
    ylimiter, xlimiter = np.max(np.abs(resid_mu)) + 0.3, [np.min(mass)-0.3, np.max(mass)+0.3]
    axs[0].set_ylim(-ylimiter, ylimiter)
    axs[1].set_ylim(-ylimiter, ylimiter)
    axs[0].set_xlim(xlimiter[0], xlimiter[1])
    axs[0].set(xlabel="Host Stellar Mass (log $M_{*}$/$[M_{\odot}]$)", ylabel='Hubble Residuals (mag)')  # Sub-plot Labels
    # axs[0].set_ylim(-1.4, 1.4); axs[0].set_xlim(7.5, 13.0)
    axs[1].get_yaxis().set_visible(False)  # Turn off y-axis labels
    axs[0].legend(loc='best')

    # Saving Figure
    if save_loc is not None:
        print('Saved figure to... ', save_loc)
        plt.savefig(save_loc)
    plt.show()
    return
def resid_v_mass_med(path, title='', save_loc=None):
    color_wheel = {'ZTF': '#D8973C', 'ATLAS': '#BD632F', 'CSP': '#72513B', 'ATLAS-ZTF': '#273E47',
                   'ZTF_SNPY': '#D8973C', 'ATLAS_SNPY': '#BD632F', 'CSP_SNPY': '#72513B', 'ATLAS-ZTF_SNPY': '#273E47',
                   'ZTF_SALT': '#D8973C', 'ATLAS_SALT': '#BD632F', 'CSP_SALT': '#72513B', 'ATLAS-ZTF_SALT': '#273E47',
                   'Pan+': 'BD632F', 'Histogram': '#3B5058',
                   '10dexfill': '#A4243B', 'mediandexfill': '#D8C99B'}
    fig, axs = plt.subplots(1, 2, figsize=(15, 6), gridspec_kw={'width_ratios': [10, 1]}, constrained_layout=True)

    # Pull data from saved text & header
    data = np.genfromtxt(path, delimiter=', ', skip_header=1, dtype=str)
    with open(path, 'r') as f:
        hdr = f.readline().split(', ')
        hdr[-1] = hdr[-1][:-1]

    # Set Arrays
    z = data[:, hdr.index('z_cmb')].astype(float)
    mass, mass_err = data[:, hdr.index('hostMass')].astype(float), data[:, hdr.index('hostMass_err')].astype(float)
    mu, mu_err = data[:, hdr.index('mu')].astype(float), data[:, hdr.index('mu_err')].astype(float)
    resid_mu, resid_mu_err = sigma_clip(mu - gen.current_cosmo().distmod(z).value, sigma=3.0), np.copy(mu_err)
    origins = data[:, hdr.index('origin')]

    # David Adjustments for Outliers
    mn,_,std = sigma_clipped_stats(mu - gen.current_cosmo().distmod(z).value,sigma=3.0)
    iGood = abs(mu - gen.current_cosmo().distmod(z).value - mn) < 3.0*std
    mu,mu_err,z,resid_mu,resid_mu_err,mass,mass_err,origins = (
        mu[iGood],mu_err[iGood],z[iGood],resid_mu[iGood],resid_mu_err[iGood],mass[iGood],mass_err[iGood],origins[iGood])

    # intrinsic dispersion added in quadrature
    mu_err = np.sqrt(mu_err ** 2.0 + 0.1 ** 2.0)

    # Make main plot
    for origin in np.unique(origins):
        format_dict = {'marker': 'o', 'fmt': 'o', 'label': origin, 'alpha': 1, 'ms': 6}
        if 'SNPY' in origin:
            format_dict['label'] = origin[:-5]
        elif 'SALT' in origin:
            format_dict['label'], format_dict['marker'] = None, '^'

        indexes = np.where(origins == origin)[0]
        axs[0].errorbar(mass[indexes], resid_mu[indexes], xerr=mass_err[indexes], yerr=resid_mu_err[indexes],
                        color=color_wheel[origin], elinewidth=0.8, **format_dict)

    # Make histogram
    axs[1].hist(resid_mu, bins=int((np.max(resid_mu) - np.min(resid_mu)) / 0.1),  # Bin Width = 0.1
                orientation="horizontal", color=color_wheel['Histogram'])

    # Extra Info
    extra_info = '$\sigma$: '+str(round(np.std(resid_mu), 4)) + ', $n$: ' + str(len(resid_mu))
    if 'merged' in path:
        extra_info += r' | SALT3: $\triangle$, SNooPy: $\bigcirc$'
    axs[0].text(np.min(mass), np.max(resid_mu),
                extra_info,
                horizontalalignment='left', verticalalignment='bottom')

    # Display Both Mass Steps
    for cut in [10, 'median']:
        num_cut, lin_color = cut, color_wheel['10dexfill']
        if cut == 'median':
            num_cut, lin_color = round(np.median(mass), 4), color_wheel['mediandexfill']

        # Get Mass Step
        mass_step_dict, resid_dict = mass_step_calc(mu, mu_err, mass, z, cut=num_cut)

        # Plot Mass Step
        lin_details = {'linestyle': '--', 'linewidth': 1.0, 'color': lin_color}
        fill_details = {'color': lin_color, 'alpha': 0.2}
        axs[0].vlines(x=num_cut, ymin=resid_dict['lower_resid']['value'], ymax=resid_dict['upper_resid']['value'],
                      **lin_details)
        axs[0].hlines(y=resid_dict['lower_resid']['value'], xmin=np.min(mass) - 0.3, xmax=num_cut,
                      **lin_details)  # Left
        axs[0].fill_between([np.min(mass) - 0.3, num_cut],
                            resid_dict['lower_resid']['value'] - resid_dict['lower_resid']['err'],
                            resid_dict['lower_resid']['value'] + resid_dict['lower_resid']['err'],
                            **fill_details)
        axs[0].hlines(y=resid_dict['upper_resid']['value'], xmin=num_cut, xmax=np.max(mass) + 0.3,
                      **lin_details)  # Right
        axs[0].fill_between([num_cut, np.max(mass) + 0.3],
                            resid_dict['upper_resid']['value'] - resid_dict['upper_resid']['err'],
                            resid_dict['upper_resid']['value'] + resid_dict['upper_resid']['err'],
                            label=str(cut) + ': ' +
                                  str(round(mass_step_dict['value'], 4)) + ' +/- ' +
                                  str(round(mass_step_dict['err'], 4)),
                            **fill_details)  # Right

    # Formatting
    fig.suptitle(title)
    axs[0].set(xlabel="Host Stellar Mass (log $M_{*}$/$[M_{\odot}]$)",
               ylabel='Hubble Residuals (mag)')  # Sub-plot Labels
    axs[1].get_yaxis().set_visible(False)  # Turn off y-axis labels
    axs[0].legend(loc='best')

    # Saving Figure
    if save_loc is not None:
        print('Saved figure to... ', save_loc)
        plt.savefig(save_loc)
    plt.show()
    return
def resid_v_mass_stacked():

    n = 0
    subplot_titles = ['Merged', 'SNooPy Only', 'SALT3 Only']
    fig, axs = plt.subplots(3, 2, figsize=(12, 6), gridspec_kw={'width_ratios': [10, 1]}, constrained_layout=True)
    for path in ['output/combiend__snpy_params_cut.txt', 'output/combiend__snpy_params_cut.txt', 'output/combiend__salt_params_cut.txt']:
        # Pull data from saved text
        data = np.genfromtxt(path, delimiter=', ', skip_header=1, dtype=str)

        # Get header
        with open(path, 'r') as f:
            hdr = f.readline().split(', ')
            hdr[-1] = hdr[-1][:-1]

        color_wheel = {'ZTF': '#81ADC8', 'ATLAS': '#EEDFAA', 'CSP': '#CD4631', 'ATLAS-ZTF': '#DEA47E',
                       'ZTF_SNPY': '#dcb160', 'ATLAS_SNPY': '#f4ac45', 'CSP_SNPY': '#af7b3f',
                       'ATLAS-ZTF_SNPY': '#694a38',
                       'ZTF_SALT': '#dcb160', 'ATLAS_SALT': '#f4ac45', 'CSP_SALT': '#af7b3f',
                       'ATLAS-ZTF_SALT': '#694a38',
                       'Pan+': 'maroon', 'Histogram': '#775A4A'}

        # Set Arrays
        z = data[:, hdr.index('z_cmb')].astype(float)
        mass, mass_err = data[:, hdr.index('hostMass')].astype(float), data[:, hdr.index('hostMass_err')].astype(float)
        mu, mu_err = data[:, hdr.index('mu')].astype(float), data[:, hdr.index('mu_err')].astype(float)
        resid_mu, resid_mu_err = sigma_clip(mu - gen.current_cosmo().distmod(z).value, sigma=3.0), np.copy(mu_err)

        # Make main plot
        for origin in np.unique(data[:, hdr.index('origin')]):
            if 'SALT' in origin:
                format_dict = {'marker': '^', 'fmt': 'o', 'label': None, 'alpha': 1, 'ms': 6}
            else:
                format_dict = {'marker': 'o', 'fmt': 'o', 'label': origin[:5], 'alpha': 1, 'ms': 5}
            indexes = np.where(data[:, hdr.index('origin')] == origin)[0]
            axs[n, 0].errorbar(mass[indexes], resid_mu[indexes], xerr=mass_err[indexes], yerr=resid_mu_err[indexes],
                            color=color_wheel[origin], **format_dict)

        # Make histogram
        axs[n, 1].hist(resid_mu, bins=10, orientation="horizontal", color=color_wheel['Histogram'])

        # Get Mass Step
        for cut in [10, 'median']:
            if cut == 'median':
                num_cut = round(np.median(mass), 4)
                lin_color, label = 'r', str(cut)
            else:
                num_cut = cut
                lin_color = 'k'

            mass_step_dict, resid_dict = mass_step_calc(mu, mu_err, mass, z, cut=num_cut)

            # Plot Mass Step
            plt_details = {'linestyle': '--', 'linewidth': 1.0, 'color': lin_color}
            fill_details = {'color': lin_color, 'alpha': 0.2}
            axs[n, 0].hlines(y=resid_dict['lower_resid']['value'], xmin=np.min(mass) - 0.3, xmax=num_cut,
                          **plt_details)  # Left
            axs[n, 0].hlines(y=resid_dict['upper_resid']['value'], xmin=num_cut, xmax=np.max(mass) + 0.3,
                          **plt_details)  # Right
            axs[n, 0].fill_between([np.min(mass) - 0.3, num_cut],
                                resid_dict['lower_resid']['value'] - resid_dict['lower_resid']['err'],
                                resid_dict['lower_resid']['value'] + resid_dict['lower_resid']['err'],
                                **fill_details)  # Left
            axs[n, 0].fill_between([num_cut, np.max(mass) + 0.3],
                                resid_dict['upper_resid']['value'] - resid_dict['upper_resid']['err'],
                                resid_dict['upper_resid']['value'] + resid_dict['upper_resid']['err'],
                                **fill_details)  # Right
            axs[n, 0].vlines(x=num_cut, ymin=resid_dict['lower_resid']['value'], ymax=resid_dict['upper_resid']['value'],
                             **plt_details)

        # Formatting
        axs[n, 0].set_title(subplot_titles[n] + '\n' +
                            'Scatter: ' + str(round(np.std(resid_mu), 2)) + ' | # of pts: ' + str(len(resid_mu))+ '\n'+str(cut))
        axs[n, 0].get_xaxis().set_visible(False);
        axs[n, 1].get_xaxis().set_visible(False)
        axs[n, 1].get_yaxis().set_visible(False)

        n += 1

    # Formatting
    axs[2, 0].get_xaxis().set_visible(True); axs[2, 1].get_xaxis().set_visible(True)
    # ylimiter, xlimiter = np.max(np.abs(resid_mu)) + 0.3, [np.min(mass) - 0.3, np.max(mass) + 0.3]
    # axs[0].set(xlabel="Host Stellar Mass (log $M_{*}$/$[M_{\odot}]$)",ylabel='Hubble Residuals (mag)')  # Sub-plot Labels
    # axs[1].get_yaxis().set_visible(False)  # Turn off y-axis labels
    axs[2, 0].legend(loc='best')

    # Saving Figure
    save_loc=None
    if save_loc is not None:
        print('Saved figure to... ', save_loc)
        plt.savefig(save_loc)
    plt.show()

    return
def norm_vs_91bg_hist(param, width=None):
    snenormIa = np.genfromtxt('txts/HiCAT_DR3_params.txt', delimiter=', ', skip_header=1, dtype=str)
    sne91bg = np.genfromtxt('output/combiend__snpy_params.txt', delimiter=', ', skip_header=1, dtype=str)
    snenormIa_hdr = ('objname, ra, dec, z, z_cmb, MJDs, MJDe, origin, mu, mu_err, st, st_err, Tmax, Tmax_err, EBVhost, '
                     'EBVhost_err, hostMass, hostMass_err').split(', ')
    sne91bg_hdr = ('objname, ra, dec, z, z_cmb, MJDs, MJDe, origin, mu, mu_err, st, st_err, Tmax, Tmax_err, EBVhost, '
                   'EBVhost_err, hostMass, hostMass_err').split(', ')

    set1 = sigma_clip(snenormIa[:, snenormIa_hdr.index(param)].astype(float), sigma=3.0)
    set2 = sigma_clip(sne91bg[:, sne91bg_hdr.index(param)].astype(float), sigma=3.0)
    if width is None:
        width = np.min([np.std(set1), np.std(set2)])

    # Plotting
    plt.figure(figsize=(12, 6))
    plt.hist(set1, color='#62BEC1', label='Normal SNe', alpha=1, bins=int((np.max(set1) - np.min(set1)) / width))
    plt.hist(set2, color='#5AD2F4', label='91bg-like', alpha=0.825, bins=int((np.max(set2) - np.min(set2)) / width))

    # Data Lines
    plt.axvline(x=np.median(set1), label='Normal SNe Median', linewidth=2.5, color='#52a1a3', linestyle=':')
    plt.axvline(x=np.median(set2), label='91bg-like Median', linewidth=2.5, color='#4bb0cc', linestyle='--')

    # Formatting
    plt.title('Normal SNeIa v. 91bg-like SNeIa -- '+param)
    plt.legend()

    # Save Plot
    # plt.savefig('saved/HiCATvPan+_'+title.split(' ')[4].lower()+'.png')
    # print('Saved to...', 'saved/HiCATvPan+_'+title.split(' ')[4].lower()+'.png')

    plt.show()
    return
def snpy_hist(hicat_params_file, norm_params_file, save_loc='', sigma = 3.0, st_width = 0.02, c_width = 0.02):
    # Open data
    hicat_data = np.genfromtxt(hicat_params_file, delimiter=', ', skip_header=1, dtype=str)
    with open(hicat_params_file, 'r') as f:
        hicat_hdr = f.readline().split(', ')
        hicat_hdr[-1] = hicat_hdr[-1][:-1]
    dr3_data = np.genfromtxt(norm_params_file, delimiter=', ', skip_header=1, dtype=str)
    with open(norm_params_file, 'r') as f:
        dr3_hdr = f.readline().split(', ')
        dr3_hdr[-1] = dr3_hdr[-1][:-1]

    fig, ax = plt.subplots(1, 2, figsize=(14, 4), constrained_layout=True)
    # ST plot
    set1 = sigma_clip(dr3_data[:, dr3_hdr.index('st')].astype(float), sigma=sigma)
    ax[0].hist(set1, int((np.max(set1) - np.min(set1)) / st_width), label='DR3', color='#5AD2F4')
    ax[0].axvline(x=np.median(set1), label=r'$\tilde{x}_{DR3}$ = '+str(round(np.median(set1),2)),
                  linewidth=2.5, color='#4bb0cc', linestyle='--')
    set2 = sigma_clip(hicat_data[:, hicat_hdr.index('st')].astype(float), sigma=sigma)
    ax[0].hist(set2, int((np.max(set2) - np.min(set2)) / st_width), label='HiCAT', color='#62BEC1', alpha=0.75)
    ax[0].axvline(x=np.median(set2), label=r'$\tilde{x}_{HiCAT}$ = '+str(round(np.median(set2),2)),
                  linewidth=2.5, color='#52a1a3', linestyle=':')
    ax[0].legend()
    ax[0].set_title('Stretch [st]')

    # EBVhost plot
    set1 = sigma_clip(dr3_data[:, dr3_hdr.index('EBVhost')].astype(float), sigma=sigma)
    ax[1].hist(set1, int((np.max(set1) - np.min(set1)) / c_width), label='DR3', color='#5AD2F4')
    ax[1].axvline(x=np.median(set1), label=r'$\tilde{x}_{DR3}$ = '+str(round(np.median(set1),2)),
                  linewidth=2.5, color='#4bb0cc', linestyle='--')
    set2 = sigma_clip(hicat_data[:, hicat_hdr.index('EBVhost')].astype(float), sigma=sigma)
    ax[1].hist(set2, int((np.max(set2) - np.min(set2)) / c_width), label='HiCAT', color='#62BEC1', alpha=0.75)
    ax[1].axvline(x=np.median(set2), label=r'$\tilde{x}_{HiCAT}$ = '+str(round(np.median(set2),2)),
                  linewidth=2.5, color='#52a1a3', linestyle=':')
    ax[1].legend()
    ax[1].set_title('Color [E(B-V) Host]')
    # ax[1].set_xlim(-0.05, 0.45)
    # ax[1].invert_xaxis()

    # plt.suptitle('HiCAT 91bg-like Type Ia SNe v. DR3 Normal Type Ia SNe\n SNooPy paramaters')
    if len(save_loc) > 0:
        plt.savefig(save_loc, dpi=300)
    plt.show()
    return
def salt_hist(hicat_params_file, norm_params_file, save_loc='', sigma = 3.0, st_width = 0.3, c_width = 0.08):
    # Open data
    hicat_data = np.genfromtxt(hicat_params_file, delimiter=', ', skip_header=1, dtype=str)
    with open(hicat_params_file, 'r') as f:
        hicat_hdr = f.readline().split(', ')
        hicat_hdr[-1] = hicat_hdr[-1][:-1]
    dr3_data = np.genfromtxt(norm_params_file, delimiter=', ', skip_header=1, dtype=str)
    with open(norm_params_file, 'r') as f:
        dr3_hdr = f.readline().split(', ')
        dr3_hdr[-1] = dr3_hdr[-1][:-1]

    fig, ax = plt.subplots(1, 2, figsize=(14, 4), constrained_layout=True)

    # x1 plot
    set1 = sigma_clip(dr3_data[:, dr3_hdr.index('x1')].astype(float), sigma=sigma)
    set2 = sigma_clip(hicat_data[:, hicat_hdr.index('x1')].astype(float), sigma=sigma)
    ax[0].hist(set1, int((np.max(set1) - np.min(set1)) / st_width), label='DR3', color='#5AD2F4')
    ax[0].hist(set2, int((np.max(set2) - np.min(set2)) / st_width), label='HiCAT', color='#62BEC1', alpha=0.75)
    ax[0].axvline(x=np.median(set1), label=r'$\tilde{x}_{DR3}$ = '+str(round(np.median(set1),2)),
                  linewidth=2.5, color='#4bb0cc', linestyle='--')
    ax[0].axvline(x=np.median(set2), label=r'$\tilde{x}_{HiCAT}$ = '+str(round(np.median(set2),2)),
                  linewidth=2.5, color='#52a1a3', linestyle=':')
    ax[0].legend()
    ax[0].set_title('Stretch [x1]')

    # c plot
    set1 = sigma_clip(dr3_data[:, dr3_hdr.index('c')].astype(float), sigma=sigma)
    set2 = sigma_clip(hicat_data[:, hicat_hdr.index('c')].astype(float), sigma=sigma)
    ax[1].hist(set1, int((np.max(set1) - np.min(set1)) / c_width), label='DR3', color='#5AD2F4')
    ax[1].hist(set2, int((np.max(set2) - np.min(set2)) / c_width), label='HiCAT', color='#62BEC1', alpha=0.75)
    ax[1].axvline(x=np.median(set1), label=r'$\tilde{x}_{DR3}$ = '+str(round(np.median(set1),2)),
                  linewidth=2.5, color='#4bb0cc', linestyle='--')
    ax[1].axvline(x=np.median(set2), label=r'$\tilde{x}_{HiCAT}$ = '+str(round(np.median(set2),2)),
                  linewidth=2.5, color='#52a1a3', linestyle=':')
    ax[1].legend()
    ax[1].set_title('Color [c]')

    # plt.suptitle('HiCAT 91bg-like Type Ia SNe v. DR3 Normal Type Ia SNe\n SALT paramaters')
    if len(save_loc) > 0:
        print('[+++] Saved figure to...', save_loc)
        plt.savefig(save_loc, dpi=300)
    plt.show()
    return
def update_readme_plots():
    # Update Mass Plots
    resid_v_mass_med(path='output/combiend__snpy_params_cut.txt',
                     title='Hubble Residual v. Host Stellar Mass of CSP-ATLAS-ZTF 91bg-like SNe Ia [SNooPy]',
                     save_loc='saved/readme_plots/csp-atlas-ztf_snpy_resid_v_mass.png')
    resid_v_mass_med(path='output/combiend__salt_params_cut.txt',
                     title='Hubble Residual v. Host Stellar Mass of CSP-ATLAS-ZTF 91bg-like SNe Ia [SALT3]',
                     save_loc='saved/readme_plots/csp-atlas-ztf_salt_resid_v_mass.png')
    resid_v_mass_med(path='output/merged_params_cut.txt',
                     title='Hubble Residual v. Host Stellar Mass of CSP-ATLAS-ZTF 91bg-like SNe Ia [SALT3-SNooPy]',
                     save_loc='saved/readme_plots/merged_resid_v_mass.png')
    resid_v_mass_med(path='txts/norm_10-11-24/normSNe_merged_10-11-24.txt',
                     title='Hubble Residual v. Host Stellar Mass of Normal SNe Ia from CSP [SNooPy]',
                     save_loc='saved/readme_plots/normIa_resid_v_mass.png')

    # Update Redshift Plots
    resid_v_z(path='output/combiend__snpy_params_cut.txt',
              title='Hubble Residual v. CMB Redshift of CSP-ATLAS-ZTF 91bg-like SNe Ia [SNooPy]',
              save_loc='saved/readme_plots/csp-atlas-ztf_snpy_resid_v_z.png')
    resid_v_z(path='output/combiend__salt_params_cut.txt',
              title='Hubble Residual v. CMB Redshift of CSP-ATLAS-ZTF 91bg-like SNe Ia [SALT3]',
              save_loc='saved/readme_plots/csp-atlas-ztf_salt_resid_v_z.png')
    resid_v_z(path='output/merged_params_cut.txt',
              title='Hubble Residual v. CMB Redshift of CSP-ATLAS-ZTF 91bg-like SNe Ia [SALT3-SNooPy]',
              save_loc='saved/readme_plots/merged_resid_v_z.png')

    # Update Histograms
    snpy_hist('output/combiend__snpy_params_cut.txt',
              'txts/norm_10-11-24/normSNe_snpy_10-11-24.txt',
              save_loc='saved/readme_plots/snpy_params_hicat_v_dr3.png')
    salt_hist('output/combiend__salt_params_cut.txt',
              'txts/norm_10-11-24/normSNe_salt_10-11-24.txt',
              save_loc='saved/readme_plots/salt_params_hicat_v_dr3.png')
    return

# Analysis Functions ------------------------------------------------------------------------------------------------ #
def mass_step_calc(mu, mu_err, mass, z, cut=10):
    if cut == 'median':
        cut = round(np.median(mass), 4)

    resid = mu - gen.current_cosmo().distmod(z).value

    upper_resid = np.average(resid[mass > cut], weights=(1/(mu_err[mass > cut]**2)))
    lower_resid = np.average(resid[mass < cut], weights=(1/(mu_err[mass < cut]**2)))

    upper_resid_err = np.std(mu_err[mass > cut]) / np.sqrt(len(mu_err[mass > cut]))  # Using Standard Error Calc
    lower_resid_err = np.std(mu_err[mass < cut]) / np.sqrt(len(mu_err[mass < cut]))

    mass_step = np.abs(upper_resid - lower_resid)
    mass_step_err = np.sqrt((lower_resid_err**2) + (upper_resid_err**2))


    return ({'value': mass_step, 'err': mass_step_err},
            {'lower_resid': {'value': lower_resid, 'err': lower_resid_err},
             'upper_resid': {'value': upper_resid, 'err': upper_resid_err}})
def dataset_analysis():
    for d_set in ['CSP', 'ATLAS', 'ZTF', 'COMBINED']:
        n_overlap = 0
        all_names, all_z, all_dec, avg_mag, avg_mag_err = [], [], [], [], []
        for algo in ['snpy', 'salt']:
            if algo == 'snpy':
                t_type = 'Tmax'
            else:
                t_type = 't0'
            sys.stdout = open(os.devnull, 'w')  # unnecessary output
            SNe = batch_load(d_set, algo)
            sys.stdout = sys.__stdout__

            for SN in SNe:
                if SN.objname not in all_names:
                    if d_set == 'COMBINED' and SN.origin == 'ATLAS-ZTF':
                        n_overlap += 1
                    all_names.append(SN.objname)
                    all_z.append(SN.z_cmb)
                    all_dec.append(SN.coords[1])
                    avg_mag.append(np.average(SN.mag))
                    time_clipped_dmag = SN.dmag[(SN.time > SN.params[t_type]['value']-5) & (SN.time < SN.params[t_type]['value']+5)]
                    if len(time_clipped_dmag) != 0:
                        avg_mag_err.append(np.average(time_clipped_dmag))

        # print('Number of SNe | Redshift | Declination | Average Mag. | Average Mag. Err (+/-5 MJD)')
        print(d_set, '&', '$'+str(len(all_names))+'$', '&',
              '$' + str(round(np.min(all_z), 4))+'$ - $'+str(round(np.max(all_z), 4))+'$', '&',
              '$' + str(round(np.min(all_dec), 6))+'$ - $'+str(round(np.max(all_dec), 6))+'$', '&',
              '$' + str(round(np.min(avg_mag), 6)) + '$ - $' + str(round(np.max(avg_mag), 6)) + '$', '&',
              '$' + str(round(np.min(avg_mag_err), 6)) + '$ - $' + str(round(np.max(avg_mag_err), 6)) + '$\\\\')
        if d_set == 'COMBINED':
            print('COMBINED Overlap:', n_overlap)
    return

# Smart Fit Functions ----------------------------------------------------------------------------------------------- #
def smart_fit_help():
    """
    Call to get examples of calls for smart_fit()
    """
    print('===========================================================================================================')
    print('Ex. Individual:', "smart_fit(fit_type='indv', data_set='CSP', algo='snpy', path='data/CSP/SN2005ke_snpy.txt')")
    print('------')
    print('Ex. Batch:', "smart_fit(fit_type='batch', data_set='CSP', algo='snpy')")
    print('------')
    print('Ex. Combined:', "smart_fit(fit_type='combiend', algo='snpy', dmag_max=1.00)")
    print('===========================================================================================================')
    return
def smart_fit(fit_type, data_set='', algo='', path=None, save_loc='', dmag_max=0.00, dflux_max=0.00):
    """
    A function to easily fit data using both algorithms.
    :param fit_type: [str] type of fitting protocol to use
    :param data_set: Default: 'CSP', [str] name of data set
    :param algo: [str] Default: 'snpy', algorithms to fit data with (snpy/salt)
    :param path: [str] Default: None, data location for individual fits (only call for fit_type='indv')
    :param save_loc: [str] location to save parameter data of SNe fits
    :param dmag_max: [float] Default: 0, magnitude error cut off for taking in data
    :param dflux_max: [float] Default: 0,  flux error cut off for taking in data
    :return: SN [list], Array of sn91bg() classes from fitting call. Returns list of 'None' if unsuccessful.
    """
    fit_type = fit_type.lower()
    if fit_type == 'indv':
        tempSN = SN91bg().make_class(data_set=data_set, path=path, dmag_max=dmag_max, dflux_max=dflux_max)
        tempSN.fit(algo)
        SNe = [tempSN]
        if SNe[0] is None:
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
        save_params_to_file(save_loc, SNe)  # Save params
        sample_cutter(save_loc, algo)  # Cut sample
    return SNe

if __name__ == '__main__':
    start = systime.time()  # Runtime tracker

    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')
