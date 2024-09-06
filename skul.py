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
from astroquery.mast import Catalogs
import scripts.general as gen
CONSTANTS = gen.get_constants()

class sn91bg():
    def __init__(self, objname=None, originalname=None, coords=(0.00, 0.00), z=0.00, origin=None):
        self.objname = objname
        self.originalname = originalname
        self.coords = coords
        self.z = z
        self.z_cmb = np.nan
        self.origin = origin

        self.period = (999999.9, 999999.9)
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
    def plot(self, y_type='mag', save_loc='', zoom=0, lines=False):
        print('[+++] Plotting LC of '+self.objname+'...')
        filter_dict = {'u': 'teal', 'g': 'green', 'r': 'red', 'i': 'indigo', 'B': 'blue',
                       'V0': 'violet', 'V1': 'purple', 'V': 'red', 'Y': 'goldenrod', 'Hdw': 'tomato', 'H': 'salmon',
                       'J': 'aquamarine', 'Jrc2': 'cadetblue', 'Jdw': 'turquoise', 'Ydw': 'olive',
                       'c': 'cyan', 'o': 'orange', 'ZTF_g': 'green', 'ZTF_r': 'red', 'ZTF_i': 'indigo'}

        # Plot
        plt.figure(figsize=(12, 6))
        for filter in np.unique(self.filters):
            indexs = np.where(self.filters == filter)[0]
            if y_type == 'mag':
                plt.errorbar(self.time[indexs], self.mag[indexs], yerr=self.dmag[indexs], fmt='o',
                             color=filter_dict[self.filters[indexs][0]], label=self.filters[indexs][0])
                if lines:
                    plt.plot(self.time[indexs], self.mag[indexs], color=filter_dict[self.filters[indexs][0]])
            elif y_type == 'flux':
                plt.errorbar(self.time[indexs], self.flux[indexs], yerr=self.dflux[indexs], fmt='o',
                             color=filter_dict[self.filters[indexs][0]], label=self.filters[indexs][0])
                if lines:
                    plt.plot(self.time[indexs], self.flux[indexs], color=filter_dict[self.filters[indexs][0]])

        # Tmax line
        if len(self.params) > 0:
            if 'Tmax' in list(self.params.keys()):
                plt.axvline(x=self.params['Tmax']['value'], color='maroon', ls='-.', label='Tmax')
            elif 't0' in list(self.params.keys()):
                plt.axvline(x=self.params['t0']['value'], color='maroon', ls='--', label='Tmax')

        # Format
        if y_type == 'mag':
            plt.gca().invert_yaxis()
        elif y_type == 'flux':
            plt.ylim(0)

        if zoom > 0:
            if 'Tmax' in list(self.params.keys()):
                plt.xlim(self.params['Tmax']['value']-zoom, self.params['Tmax']['value']+zoom)
            elif 't0' in list(self.params.keys()):
                plt.xlim(self.params['t0']['value']-zoom, self.params['t0']['value']+zoom)

        plt.title('Lightcurve -- '+self.objname+' | '+self.originalname+' -- '+y_type)
        plt.xlabel('MJD');
        plt.ylabel(y_type)
        plt.legend()
        if len(save_loc) > 0:
            plt.savefig(save_loc + obj + '_lc.png')
            print(self.objname, '-- Plot saved to', save_loc + self.objname + '_lc.png')
        plt.show()
        systime.sleep(2)

        return
    def subplots(self, y_type='mag', save_loc='', zoom=0, ticks=5, lines=False):
        print('[+++] Plotting LC of '+self.objname+'...')
        filter_dict = {'u': 'teal', 'g': 'green', 'r': 'red', 'i': 'indigo', 'B': 'blue',
                       'V0': 'violet', 'V1': 'purple', 'V': 'red', 'Y': 'goldenrod', 'Hdw': 'tomato', 'H': 'salmon',
                       'J': 'aquamarine', 'Jrc2': 'cadetblue', 'Jdw': 'turquoise', 'Ydw': 'olive',
                       'c': 'cyan', 'o': 'orange', 'ZTF_g': 'green', 'ZTF_r': 'red', 'ZTF_i': 'indigo'}

        # Plot
        unique_filters, num_plts = np.unique(self.filters), len(np.unique(self.filters))
        if self.origin == 'CSP':
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
            plt.subplot(size[0], size[1], i+1)
            indexs = np.where(self.filters == unique_filters[i])[0]
            if y_type == 'mag':
                plt.errorbar(self.time[indexs], self.mag[indexs], yerr=self.dmag[indexs], fmt='o', ms=4, elinewidth=0.3,
                             color=filter_dict[self.filters[indexs][0]], label=self.filters[indexs][0])
            elif y_type == 'flux':
                plt.errorbar(self.time[indexs], self.flux[indexs], yerr=self.dflux[indexs], fmt='o', ms=4, elinewidth=0.3,
                             color=filter_dict[self.filters[indexs][0]], label=self.filters[indexs][0])

            # Tmax line
            if len(self.params) > 0:
                if 'Tmax' in list(self.params.keys()):
                    plt.axvline(x=self.params['Tmax']['value'], color='maroon', ls='-.', label='Tmax')
                elif 't0' in list(self.params.keys()):
                    plt.axvline(x=self.params['t0']['value'], color='maroon', ls='-.', label='Tmax')

            # Format
            if y_type == 'mag':
                plt.gca().invert_yaxis()
            elif y_type == 'flux':
                plt.ylim(0)
            if i > 0:
                plt.gca().get_yaxis().set_visible(False) # Turn off y-axis labels
            # plt.xticks(np.arange(np.min(self.time), np.max(self.time), ((np.max(self.time) - np.min(self.time)) / ticks)))

            if zoom > 0:
                if 'Tmax' in list(self.params.keys()):
                    plt.xlim(self.params['Tmax']['value']-zoom, self.params['Tmax']['value']+zoom)
                elif 't0' in list(self.params.keys()):
                    plt.xlim(self.params['t0']['value']-zoom, self.params['t0']['value']+zoom)

            plt.xlabel('MJD'); plt.legend() # plt.ylabel(y_type);
        plt.suptitle('Lightcurve -- '+self.objname+' | '+self.originalname+' -- '+y_type)
        if len(save_loc) > 0:
            plt.savefig(save_loc + obj + '_lc.png')
            print(self.objname, '-- Plot saved to', save_loc + self.objname + '_lc.png')
        plt.show()
        systime.sleep(2)
        return
    def save_class(self, save_loc):
        print('[+++] '+self.objname+' -- Saving class to '+save_loc+'classes/'+self.objname+'_class.txt')
        with open(save_loc+'classes/'+self.objname+'_class.txt', 'w') as f:
            f.write(self.origin + ' ' + self.objname + ' ' + self.originalname +
                    ' ' + str(self.coords[0]) + ' ' + str(self.coords[1]) +
                    ' ' + str(self.z) + ' ' + str(self.z_cmb) + '\n')
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
            self.origin, self.objname, self.originalname, self.coords, self.z, self.z_cmb = (
                hdr[0], hdr[1], hdr[2], (float(hdr[3]), float(hdr[4])), float(hdr[5]), float(hdr[6][:-1]))

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
        return
    def set_data(self, zp, filters, time, flux, dflux, mag, dmag):
        self.zp = np.copy(zp)
        self.filters = np.copy(filters)
        self.time = np.copy(time)
        self.flux = np.copy(flux)
        self.dflux = np.copy(dflux)
        self.mag = np.copy(mag)
        self.dmag = np.copy(dmag)
        return
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
    def snpy_fit(self, save_loc, show_plot=True, quiet=False): # use_saved=False, snpy_plots=True, save_plots=True,
        print('[+++] '+self.objname+' -- Fitting data with SNooPy...')
        load_path = save_loc + 'ascii/' + self.objname + '_snpy.txt'
        save_path = save_loc + 'models/' + self.objname + '_EBV_model2.snpy'

        # Check quiet
        if quiet:
            sys.stdout = open(os.devnull, 'w')

        # Load Data
        try:
            n_s = snpy.get_sn(load_path)
        except:
            self.params.update({'mu': {'value': 0.00, 'err': 0.00}})
            print('[!!!] Failed to load ASCII file')
            return
        n_s.k_version = '91bg'
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

                param_names = ['mu', 'st', 'Tmax', 'EBVhost']
                snpy_param_names = ['DM', 'st', 'Tmax', 'EBVhost']
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
    def salt_fit(self, save_loc, show_plot=True):
        print('[+++] '+self.objname+' -- Fitting data with SALT3...')
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
            model.set(z=self.z, t0=self.time[np.where(self.flux == np.max(self.flux))[0][0]])  # set the model's redshift.
            result, fitted_model = sncosmo.fit_lc(data, model, ['t0', 'x0', 'x1', 'c'], bounds={'x1': (-5, 5)})

            param_names = ['t0', 'x0', 'x1', 'c']
            for i in range(len(param_names)):
                self.params.update({param_names[i]: {'value': result.parameters[i+1],
                                                     'err': result.errors[param_names[i]]}})

            # Calculate
            pho_mB = -2.5 * np.log10(self.params['x0']['value']) + mB_const
            pho_mB_err = np.abs((-2.5 * self.params['x0']['err']) / (self.params['x0']['err'] * np.log(10)))

            mu = pho_mB + (alpha * self.params['x1']['value']) - (beta * self.params['c']['value']) - M0
            mu_err = np.sqrt(pho_mB_err ** 2 + (np.abs(alpha) * self.params['x1']['err']) ** 2 + (
                        np.abs(alpha) * self.params['c']['err']) ** 2)

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
                self.params.update({'hostMass': {'value': mass_key[self.objname]['mass'],
                                                 'err': mass_key[self.objname]['mass_err']}})
                print('[+++] Mass taken from mass key!')
                return

        # Getting host data -- checks GLADE then PANSTARRS
        gMag = -999.00
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
                    print(host_data['l'].loc[0], host_data['b'].loc[0])
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
        print('      Success!', self.objname, 'host galaxy has a mass of:', host_mass, '+/-', host_mass_err, 'logM_* / [Msun]')

        # Update mass key
        if use_key:
            with open(CONSTANTS['mass_key_txt'], 'a') as f:
                print('      Updating mass key with ' + self.objname + '...')
                f.write(self.objname + ', ' + str(host_mass) + ', ' + str(host_mass_err) + '\n')

        print('      Removing GHOST data...')
        shutil.rmtree('default/ghost_stuff/')  # Clear messy data
        os.mkdir('default/ghost_stuff/')
        return

def save_params_to_file(save_loc, SNe):
    print('[+++] Saving params to '+save_loc+'...')
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
def class_creation(data_set, path, dmag_max=0, dflux_max=0):
    if data_set == 'CSP':
        with open(path, 'r') as f:
            objname, z, ra, dec = f.readline().split(' ')
            originalname = path.split('/')[-1].split('_')[0]
        objname, ra, dec, z = objname[2:], float(ra), float(dec[:-1]), float(z)
        tempSN = sn91bg(objname, originalname, (ra, dec), z, 'CSP')

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

                    # Calculate Magnitude
                    # n_flux = 10 ** ((n_mag - zp) / -2.5)
                    # n_dflux = 10 ** (n_mag / -2.5)

                    n_flux = 10 ** ((n_mag - zp) / -2.5)
                    n_dflux = np.abs(n_flux)*np.log(10)*((1/2.5)*n_dmag)

                    tempSN.zp = np.append(tempSN.zp, zp)
                    tempSN.filters = np.append(tempSN.filters, filter)
                    tempSN.time = np.append(tempSN.time, n_time)
                    tempSN.flux = np.append(tempSN.flux, n_flux)
                    tempSN.dflux = np.append(tempSN.dflux, n_dflux)
                    tempSN.mag = np.append(tempSN.mag, n_mag)
                    tempSN.dmag = np.append(tempSN.dmag, n_dmag)

        tempSN.clean_data(dmag_max, dflux_max)
    elif data_set == 'ATLAS':
        data = np.genfromtxt(path, delimiter=',', dtype=str, skip_header=1)
        if len(data) == 0:
            print('[!!!] File [' + path + '] empty!')
            return None
        ra, dec = np.average(data[:, 1].astype(float)), np.average(data[:, 2].astype(float))
        objname, z = gen.TNS_objname_z(ra, dec)
        originalname = path.split('/')[-1].split('.')[0]
        z = np.nan if z == 'None' else float(z)
        tempSN = sn91bg(objname, originalname, (ra, dec), z, 'ATLAS')
        tempSN.set_data(zp=data[:, 7].astype(float), filters=data[:, 6], time=data[:, 8],
                        flux=data[:, 16], dflux=data[:, 17], mag=data[:, 3], dmag=data[:, 4])
        tempSN.clean_data(dmag_max, dflux_max)
    elif data_set == 'ZTF':
        data = np.genfromtxt(path, delimiter=None, dtype=str, skip_header=56)
        if len(data) == 0:
            print('[!!!] File [' + path + '] empty!')
            return None
        with (open(path, 'r') as f):
            # ztf_spread = 200
            ztf_spread = float(CONSTANTS['ztf_spread'])

            hdr = f.readlines()
            ra, dec = float(hdr[3].split(' ')[-2]), float(hdr[4].split(' ')[-2])
            objname, z = gen.TNS_objname_z(ra, dec)
            originalname = path.split('/')[-1].split('.')[0].split('_')[1]
            z = np.nan if z == 'None' else float(z)

            # Get magnitudes m = -2.5log(F) + zp
            time, flux, dflux, zp, filters = data[:, 22], data[:, 24], data[:, 25], data[:, 20], data[:, 4]
            valid_ints = np.unique(np.hstack((np.where(flux != 'null')[0], np.where(dflux != 'null')[0])))
            time, zp, filters = time[valid_ints].astype(float), zp[valid_ints].astype(float), filters[valid_ints]
            flux, dflux = flux[valid_ints].astype(float), dflux[valid_ints].astype(float)
            time = time - 2400000.5 # JD to MJD
            mag, dmag = ((-2.5 * np.log10(flux)) + zp), (2.5 * np.log10(dflux))

            # Adjusting around tmax
            if ztf_spread != 0 and len(time) != 0:
                t_max_guess = time[np.where(flux == np.max(flux))[0][0]]
                zoom_indexes = np.where(time[np.where(time > t_max_guess - ztf_spread)[0]] < t_max_guess + ztf_spread)[0]
                time = time[zoom_indexes]
                flux = flux[zoom_indexes]
                dflux = dflux[zoom_indexes]
                mag = mag[zoom_indexes]
                dmag = dmag[zoom_indexes]
                zp = zp[zoom_indexes]
                filters = filters[zoom_indexes]

        if len(time) == 0:
            print('[!!!] File has no valid values!')
            return None

        tempSN = sn91bg(objname, originalname, (ra, dec), z, 'ZTF')
        tempSN.set_data(zp=zp, filters=filters, time=time,
                        flux=flux, dflux=dflux, mag=mag, dmag=dmag)
        tempSN.clean_data(dmag_max, dflux_max)
    else:
        raise ValueError("Data set '" + data_set + "' not recognized")
    return tempSN
# ------------------------------------------------------------------------------------------------------------------- #
def combined_fit(algo='snpy', cut=False, dmag_max=0, dflux_max=0):
    sys.stdout = open(os.devnull,'w') # Lots of unecessary output
    csp_files = glob.glob('data/CSP/*.txt')
    CSP_SNe = {}
    for file in csp_files:
        tempSN = class_creation('CSP', file, dmag_max, dflux_max)
        if tempSN is not None:
            CSP_SNe.update({tempSN.objname: tempSN})
    atlas_files = glob.glob('data/ATLAS/*.txt')
    ATLAS_SNe, atlas_names = {}, []
    for file in atlas_files:
        tempSN = class_creation('ATLAS', file, dmag_max, dflux_max)
        if tempSN is not None:
            ATLAS_SNe.update({tempSN.objname: tempSN})
    ztf_files = glob.glob('data/ZTF/*.txt')
    ZTF_SNe = {}
    for file in ztf_files:
        tempSN = class_creation('ZTF', file, dmag_max, dflux_max)
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
            ATLAS_SNe[name].set_data(new_zp, new_filters, new_time, new_flux, new_dflux, new_mag, new_dmag)
            combined_SNe.append(ATLAS_SNe[name])
        elif name in list(CSP_SNe.keys()):
            combined_SNe.append(CSP_SNe[name])
        elif name in list(ATLAS_SNe.keys()):
            combined_SNe.append(ATLAS_SNe[name])
        elif name in list(ZTF_SNe.keys()):
            combined_SNe.append(ZTF_SNe[name])

    # Fitting data
    fit_combined_SNe = []
    if algo == 'snpy':
        for SN in combined_SNe:
            print('[', combined_SNe.index(SN)+1, '/', len(combined_SNe), '] Fitting data with SNooPy for '+SN.objname+'...')
            print('-----------------------------------------------------------------------------------------------')
            SN.write_snpy_ASCII(save_loc=CONSTANTS['combined_saved_loc'] + 'ascii/')
            SN.snpy_fit(save_loc=CONSTANTS['combined_saved_loc'])
            if SN.params['mu']['value'] <= 0.00:
                continue
            SN.get_host_mass()
            SN.save_class(CONSTANTS['combined_saved_loc'])
            if 'hostMass' in list(SN.params.keys()):
                fit_combined_SNe.append(SN)
    elif algo == 'salt':
        for SN in combined_SNe:
            print('[', combined_SNe.index(SN)+1, '/', len(combined_SNe), '] Fitting data with SALT3 for '+SN.objname+'...')
            print('-----------------------------------------------------------------------------------------------')
            SN.salt_fit(save_loc=CONSTANTS['salt_combined_loc'])
            if SN.params['mu']['value'] <= 0.00:
                return None
            SN.get_host_mass()
            SN.save_class(CONSTANTS['salt_combined_loc'])
            if 'hostMass' in list(tempSN.params.keys()):
                fit_combined_SNe.append(SN)
    else:
        raise ValueError("[!!!] Invalid algorithm selected ['snpy'/'salt']")
    print('Sucessfully fit [', len(fit_combined_SNe), '/', len(combined_SNe), ']!')

    if algo == 'snpy':
        save_params_to_file(CONSTANTS['combined_saved_loc']+'combined_params.txt', fit_combined_SNe)
    elif algo == 'salt':
        save_params_to_file(CONSTANTS['salt_combined_loc']+'combined_params.txt', fit_combined_SNe)

    return fit_combined_SNe
def batch_fit(data_set, algo='snpy', cut=False, dmag_max=0, dflux_max=0):
    SNe, files = [], glob.glob('data/' + data_set + '/*.txt')
    for path in files:
        print('[', files.index(path) + 1, '/', len(files), ']')
        print('-----------------------------------------------------------------------------------------------')
        tempSN = indvisual_fit(data_set, path, algo, dmag_max, dflux_max)
        if tempSN is not None and 'hostMass' in list(tempSN.params.keys()):
            SNe.append(tempSN)
    print('Sucessfully fit [', len(SNe), '/', len(files), ']!')

    # Cut sample
    if cut:
        SNe = sample_cutter(SNe, data_set, algo=algo)

    # Save params to file
    if algo == 'snpy':
        save_params_to_file(CONSTANTS[data_set.lower() + '_saved_loc'] + data_set.lower() + '_params.txt',
                            SNe)
    if algo == 'salt':
        save_params_to_file(CONSTANTS['salt_'+data_set.lower()+'_loc'] + data_set.lower() + '_params.txt',
                            SNe)
    return SNe
def indvisual_fit(data_set, path, algo='snpy', dmag_max=0, dflux_max=0):
    print(path)
    if algo == 'snpy':
        tempSN = class_creation(data_set, path, dmag_max, dflux_max)
        if tempSN is None:
            return None
        tempSN.write_snpy_ASCII(save_loc=CONSTANTS[data_set.lower() + '_saved_loc'] + 'ascii/')
        tempSN.snpy_fit(save_loc=CONSTANTS[data_set.lower() + '_saved_loc'])
        if tempSN.params['mu']['value'] <= 0.00:
            return None
        tempSN.get_host_mass(use_key=True)
        tempSN.save_class(CONSTANTS[data_set.lower() + '_saved_loc'])
    elif algo == 'salt':
        tempSN = class_creation(data_set, path)
        if tempSN is None:
            return None
        tempSN.salt_fit(save_loc=CONSTANTS['salt_' + data_set.lower() + '_loc'])
        if tempSN.params['mu']['value'] <= 0.00:
            return None
        tempSN.get_host_mass()
        if tempSN.params['mu']['value'] <= 0.00:
            return None
        tempSN.save_class(CONSTANTS['salt_' + data_set.lower() + '_loc'])
    else:
        raise ValueError("[!!!] Invalid algorithm selected ['snpy'/'salt']")
    return tempSN
# ------------------------------------------------------------------------------------------------------------------- #
def indivisual_load(path):
    tempSN = sn91bg()
    tempSN.load_from_file(path)
    return tempSN
def batch_load(data_set, algo='snpy'):
    SNe = []
    for path in glob.glob('saved/' + algo + '/' + data_set.lower() + '/classes/*_class.txt'):
        SNe.append(indivisual_load(path))
    return SNe
# ------------------------------------------------------------------------------------------------------------------- #
def sample_cutter(SNe, data_set, algo='snpy', save=True):
    print('[+++] Cutting sample...')
    new_SNe = []
    if algo == 'snpy':
        print('[+++] Cutting sample for SNooPy data...')
        for SN in SNe:
            readout = '[' + str(SNe.index(SN) + 1) + '/' + str(len(SNe)) + '] -- ' + SN.objname + ' '
            resid = float(SN.params['mu']['value']) - gen.current_cosmo().distmod(float(SN.z)).value
            resid -= np.median(resid)

            if float(SN.params['EBVhost']['value']) < -1 or float(SN.params['EBVhost']['value']) > 1:
                readout += '-- EBVhost out of range! -- ' + str(SN.params['EBVhost']['value'])
                print(readout)
                continue
            if float(SN.params['EBVhost']['err']) > 0.5:
                readout += '-- EBVhost errors out of range! -- ' + str(SN.params['EBVhost']['err'])
                print(readout)
                continue

            if float(SN.params['st']['value']) < 0.3 or float(SN.params['st']['value']) > 1.0:
                readout += '-- Stretch out of range! -- ' + str(SN.params['st']['value'])
                print(readout)
                continue
            if float(SN.params['st']['err']) > 0.1:
                readout += '-- Stretch error out of range! -- ' + str(SN.params['st']['err'])
                print(readout)
                continue

            if float(SN.params['Tmax']['err']) > 1:
                readout += '-- Maximum time error out of range! -- ' + str(SN.params['Tmax']['err'])
                print(readout)
                continue

            if float(SN.z_cmb) < 0.015:
                readout += '-- Redshift out of range! -- ' + str(SN.z_cmb)
                print(readout)
                continue

            if float(SN.params['hostMass']['value']) <= 0.00:
                readout += '-- Host mass null! -- ' + str(SN.params['hostMass']['value'])
                print(readout)
                continue

            new_SNe.append(SN)

        print('Successfully cut data [', len(new_SNe), '/', len(SNe), '] !')
    elif algo == 'salt':
        print('[+++] Cutting SALT3 results...')
        for SN in SNe:
            readout = '[' + str(SNe.index(SN) + 1) + '/' + str(len(SNe)) + '] -- ' + SN.objname + ' '

            if float(SN.params['x1']['value']) < -3 or float(SN.params['x1']['value']) > 3:
                readout += '-- X1 out of range!'
                print(readout)
                continue
            if float(SN.params['x1']['err']) > 1:
                readout += '-- X1 error too large!'
                print(readout)
                continue

            if float(SN.params['c']['value']) < -0.3 or float(SN.params['c']['value']) > 0.3:
                readout += '-- c out of range!'
                print(readout)
                continue
            if float(SN.params['c']['err']) > 0.1:
                readout += '-- c error too large!'
                print(readout)
                continue

            if float(SN.params['t0']['err']) > 2:
                readout += '-- t0 error too large!'
                print(readout)
                continue

            new_SNe.append(SN)

        print('Successfully cut data [', len(new_SNe), '/', len(SNe), '] !')

    # Save params to file
    if save:
        if algo == 'snpy':
            save_params_to_file(CONSTANTS[data_set.lower() + '_saved_loc'] + data_set.lower() + '_cut_params.txt',
                                SNe)
        elif algo == 'salt':
            save_params_to_file(CONSTANTS['salt_'+data_set.lower()+'_loc'] + data_set.lower() + '_cut_params.txt',
                                SNe)

    return new_SNe
def file_sample_cutter(path, algo, save=True):
    print('[+++] Cutting sample from file...')

    valid_lines = []
    chi2 = []
    with open(path, 'r') as f:
        hdr = f.readline().split(', ')
        hdr[-1] = hdr[-1][:-1]
        data = f.readlines()
        uncut_num = len(data)
        for line in data:
            n_line = line.split(', ')
            if algo == 'snpy':
                if float(n_line[hdr.index('EBVhost')]) < -0.2 or float(n_line[hdr.index('EBVhost')]) > 0.2:
                    continue
                if float(n_line[hdr.index('EBVhost_err')]) > 0.1:
                    continue
                if float(n_line[hdr.index('st')]) < 0.3 or float(n_line[hdr.index('EBVhost')]) > 1.0:
                    continue
                if float(n_line[hdr.index('st_err')]) > 0.1:
                    continue
                if float(n_line[hdr.index('Tmax_err')]) > 1:
                    continue
                if float(n_line[hdr.index('chi2_err')]) > 2:
                    continue
            elif algo == 'salt':
                if float(n_line[hdr.index('x1')]) < -3 or float(n_line[hdr.index('x1')]) > 3:
                    continue
                if float(n_line[hdr.index('x1_err')]) > 1:
                    continue
                if float(n_line[hdr.index('c')]) < -0.3 or float(n_line[hdr.index('c')]) > 0.3:
                    continue
                if float(n_line[hdr.index('c_err')]) > 0.1:
                    continue
                if float(n_line[hdr.index('t0')]) > 2:
                    continue

            if float(n_line[hdr.index('z')]) < 0.015:
                continue
            if float(n_line[hdr.index('hostMass')]) <= 0:
                continue
            valid_lines.append(line)
    print('Successfully cut data [', len(valid_lines), '/', uncut_num, '] !')

    # Save params to file
    if save:
        print('Saving to... ', path[:-16]+path[-14:])
        with open(path[:-16]+path[-14:], 'w') as f:
            f_hdr = hdr[0]
            for n in range(1, len(hdr)):
                f_hdr += ', ' + hdr[n]
            f_hdr += '\n'
            f.write(f_hdr)
            for line in valid_lines:
                f.write(line)
    return
# ------------------------------------------------------------------------------------------------------------------- #
def residual_plotter(path, x_params='Redshift', labels=False):
    # Pull data from saved text
    data = np.genfromtxt(path, delimiter=', ', skip_header=1, dtype=str)
    data_set = path.split('/')[-1].split('_')[0].upper()
    algo = path.split('_')[-3].upper()

    fig, axs = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'width_ratios': [10, 1]}, constrained_layout=True)
    color_wheel = {'ZTF': '#81ADC8', 'ATLAS': '#EEDFAA', 'CSP': '#CD4631', 'ATLAS-ZTF': '#DEA47E'}
    key = {'ZTF': True, 'ATLAS': True, 'CSP': True, 'ATLAS-ZTF': True}

    # Plot points
    for origin in np.unique(data[:, 7]):
        indexs = np.where(data[:, 7] == origin)[0]
        resid_mu = data[:, 8].astype(float)[indexs] - gen.current_cosmo().distmod(data[:, 4].astype(float)[indexs]).value
        resid_mu_err = data[:, 9].astype(float)[indexs]

        if x_params == 'Host Mass':
            x_axis, x_axis_err = data[:, 15].astype(float)[indexs], data[:, 16].astype(float)[indexs]
        elif x_params == 'Redshift':
            x_axis, x_axis_err = data[:, 4].astype(float)[indexs], None
        else:
            raise ValueError("[!!!] Invalid x_params ['Host Mass'/'Redshift']")
        axs[0].errorbar(x_axis, resid_mu, xerr=x_axis_err, yerr=resid_mu_err,
                        color=color_wheel[origin], fmt='o', label=origin)
        # Labels
        if labels:
            for i in range(len(resid_mu)):
                axs[0].text(x_axis[i], resid_mu[i], str(data[:, 0][indexs][i]), size='x-small', va='top')

    # Plot histogram
    all_resid_mu = data[:, 8].astype(float) - gen.current_cosmo().distmod(data[:, 4].astype(float)).value
    axs[1].hist(all_resid_mu, bins=40, orientation="horizontal", color='#9E6240')

    # Formatting
    ylimiter = np.max(np.abs(resid_mu))+0.5
    axs[0].set_ylim(-ylimiter, ylimiter); axs[1].set_ylim(-ylimiter, ylimiter)
    fig.suptitle("Hubble Residuals vs. " + x_params + " of '"+data_set+"' 91bg-like SNe Ia -- "+algo, fontsize=15)

    fig.suptitle("Hubble Residuals vs. " + x_params + " of '"+data_set+"' 91bg-like SNe Ia\n" +  # Figure Title
             ' | Scatter: ' + str(round(np.std(all_resid_mu), 2)) + ' | # of pts: ' + str(len(all_resid_mu)), size='medium')

    axs[0].set(xlabel=x_params, ylabel='Hubble Residuals')  # Sub-plot Labels
    axs[1].get_yaxis().set_visible(False) # Turn off y-axis labels
    axs[0].legend()
    print('Saved figure to... ', path[:-len(path.split('/')[-1])]+data_set.lower()+'_resid_v_'+x_params.lower()+'.png')
    plt.savefig(path[:-len(path.split('/')[-1])]+data_set.lower()+'_resid_v_'+x_params.lower()+'.png')
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
    print('Saved figure to... ', save_loc+path.split('/')[-1].split('_')[0]+'_hist.png')
    plt.savefig(save_loc+path.split('/')[-1].split('_')[0]+'_hist.png')
    plt.show()
    return
# ------------------------------------------------------------------------------------------------------------------- #
def all_fit(plot=True):
    for algo_sel in ['snpy', 'salt']:
        combined_fit(algo='salt', cut=False)
        for source in ['CSP', 'ATLAS', 'ZTF']:
            batch_fit(source, algo=algo_sel, cut=True)

    if plot:
        for algo_sel in ['snpy', 'salt']:
            for source in ['csp', 'atlas', 'ztf', 'combined']:
                n_path = 'saved/'+algo_sel+'/'+source+'/'+source+'_params.txt'
                residual_plotter(path=n_path, x_params='Redshift', labels=False)
                residual_plotter(path=n_path, x_params='Host Mass', labels=False)
                histogram_plotter(path=n_path, param_bins=[5, 5, 5, 5, 5])
    return
def residual_v_EBVhost(data_loc, val_cut = 0.2, err_cut = 0.1):
    data = np.genfromtxt(data_loc, dtype=str, skip_header=1, delimiter=', ')
    names, origins = data[:, 0], data[:, 6]
    redshifts, mus, mu_errs = data[:, 3].astype(float), data[:, 7].astype(float), data[:, 8].astype(float)
    EBVhosts, EBVhost_errs = data[:, 13].astype(float), data[:, 14].astype(float)
    color_wheel = {'ZTF': '#81ADC8', 'ATLAS': '#EEDFAA', 'CSP': '#CD4631', 'ATLAS-ZTF': '#DEA47E'}

    plt.figure(figsize=(10, 6))
    resid_arr = np.array([])
    for i in range(len(mus)):
        name = names[i]
        origin = origins[i]
        n_resid = mus[i] - CURRENT_COSMO.distmod(redshifts[i]).value
        n_resid_err = mu_errs[i]
        n_EBVhost = EBVhosts[i]
        n_EBVhost_err = EBVhost_errs[i]

        if n_EBVhost_err > 5:
            continue

        if val_cut != 0 and err_cut != 0:
            if n_EBVhost > val_cut or n_EBVhost < -val_cut or n_EBVhost_err > err_cut:
                continue

        resid_arr = np.append(resid_arr, n_resid)
        plt.errorbar(n_EBVhost, n_resid, xerr=n_EBVhost_err, yerr=n_resid_err, fmt='o',
                     elinewidth=1, color=color_wheel[origin])
        # plt.text(n_EBVhost, n_resid, name, size='x-small', va='top')

    plt.xlabel('EBVhost'); plt.ylabel('Residual')
    title = 'EBVhost v. Residual \n'
    title += data_loc.split('/')[-1].split('_')[-2]
    title += ' | # = ' + str(len(resid_arr)) + ' | Residual Scater = ' + str(round(np.std(resid_arr), 4))
    if val_cut != 0 and err_cut != 0:
        title += '\n' + str(val_cut) + ' < EBVhost < ' + str(val_cut) + ' | EBVhosterr < ' + str(err_cut)
    # plt.xlim(-0.2, 0.2); plt.ylim(-1, 1)
    plt.title(title)
    plt.show()
    return
def snpy_v_salt_mu():
    for d_set in ['csp', 'atlas', 'ztf']:
        snpy_data = np.genfromtxt('../output/batch_'+d_set+'_snpy_uncut_params.txt', skip_header=1, dtype=str, delimiter=', ')
        snpy_names, snpy_mu = snpy_data[:, 0], snpy_data[:, 7].astype(float)
        salt_data = np.genfromtxt('../output/batch_'+d_set+'_salt_uncut_params.txt', skip_header=1, dtype=str, delimiter=', ')
        salt_names, salt_mu = salt_data[:, 0], salt_data[:, 15].astype(float)

        all_mu = {}
        all_salt_mu, all_snpy_mu, resid = np.array([]), np.array([]), np.array([])
        for name in snpy_names:
            if name in salt_names:
                all_mu.update({name: {'salt_mu': salt_mu[np.where(salt_names == name)[0][0]],
                                      'snpy_mu': snpy_mu[np.where(snpy_names == name)[0][0]]}})
                all_snpy_mu = np.append(all_snpy_mu, snpy_mu[np.where(snpy_names == name)[0][0]])
                all_salt_mu = np.append(all_salt_mu, salt_mu[np.where(salt_names == name)[0][0]])
                resid = np.append(resid, snpy_mu[np.where(snpy_names == name)[0][0]] - salt_mu[np.where(salt_names == name)[0][0]])
        plt.title(d_set)
        plt.hist(resid)
        plt.xlim(-max(abs(resid)), max(abs(resid)))
        plt.show()
    return
def param_compare():
    dr3_data = np.genfromtxt('../txts/DR3_fits.dat', dtype=str, skip_header=1)
    data = np.genfromtxt('../output/COMBINED_combined_snpy_uncut_rv25_params.txt', dtype=str, skip_header=1, delimiter=', ')

    overlapped_st_dr3, overlapped_ebv_dr3 = np.array([]), np.array([])
    overlapped_st, overlapped_ebv = np.array([]), np.array([])
    for i in range(len(data[:, 0])):
        if 'SN'+data[i, 0] in dr3_data[:, 0]:
            overlapped_st_dr3 = np.append(overlapped_st_dr3, float(dr3_data[np.where(dr3_data[:, 0] == 'SN'+data[i, 0])[0][0], 2]))
            overlapped_ebv_dr3 = np.append(overlapped_ebv_dr3, float(dr3_data[np.where(dr3_data[:, 0] == 'SN'+data[i, 0])[0][0], 25]))
            overlapped_st = np.append(overlapped_st, float(data[i, 9]))
            overlapped_ebv = np.append(overlapped_ebv, float(data[i, 13]))
            # print('SN'+data[i, 0], dr3_data[np.where(dr3_data[:, 0] == 'SN'+data[i, 0])[0][0], 0])

    # # ST Hist
    # plt.figure(figsize=(12, 7))
    # plt.hist(overlapped_st, bins=10, color='r', label='SNPY')
    # plt.hist(overlapped_st_dr3, bins=10, color='b', label='DR3')
    # plt.title('SNPY v. DR3 Params -- st')
    # plt.legend()
    # plt.show()

    # EBVhost Hist
    plt.figure(figsize=(12, 7))
    plt.hist(overlapped_ebv, bins=30, color='r', label='SNPY')
    plt.hist(overlapped_ebv_dr3, bins=30, color='b', label='DR3')
    plt.title('SNPY v. DR3 Params -- EBVhost')
    plt.legend()
    plt.show()

    return
def param_hist_compare(set1, set2, bin_width, xrange=None, title=''):
    plt.figure(figsize=(9, 4))
    plt.title(title)
    plt.hist(set1, color='#62BEC1', label=title.split(' ')[0], alpha=0.75,
             bins=int((np.max(set1) - np.min(set1)) / bin_width))
    plt.hist(set2, color='#5AD2F4', label=title.split(' ')[2], alpha=0.75,
             bins=int((np.max(set2) - np.min(set2)) / bin_width))
    plt.xlim(xrange)
    plt.legend()
    plt.show()
    return
# ------------------------------------------------------------------------------------------------------------------- #
def output(fit_type, data_set, algo, cut=False, path=None, special_text=''):
    if fit_type == 'indv':
        SNe = [indvisual_fit(data_set, path, algo=algo)]
        if SNe[0] == None:
            print('[!!!] Fit failed!')
            return [None]
    elif fit_type == 'batch':
        SNe = batch_fit(data_set, algo=algo, dmag_max=1)
    elif fit_type == 'COMBINED':
        SNe = combined_fit(algo=algo)
    else:
        raise ValueError("Invalid type selected ['indv'/'batch'/'COMBINED']")
    # Saving
    save_params_to_file('output/'+fit_type+'_'+data_set.lower()+'_'+algo+'_uncut_'+special_text+'params.txt', SNe)
    if cut:
        SNe = sample_cutter(SNe, data_set, algo=algo, save=False)
        if len(SNe) > 0:
            save_params_to_file('output/'+fit_type+'_'+data_set.lower()+'_'+algo+'_cut_'+special_text+'params.txt', SNe)
    return SNe

if __name__ == '__main__':
    start = systime.time() # Runtime tracker

    # output('COMBINED', 'COMBINED', 'snpy', cut=True) # path='data/ATLAS/1041641640121059200.txt'

    residual_plotter('output/COMBINED_combined_snpy_uncut_params.txt', x_params='Redshift', labels=False)
    residual_plotter('output/COMBINED_combined_snpy_uncut_params.txt', x_params='Host Mass', labels=False)
    # histogram_plotter(path='saved/salt/atlas/atlas_params.txt', param_bins=[5, 5, 5, 5, 5]) # path='saved/snpy/csp/csp_params.txt', param_bins=[5, 5, 5, 5, 5]

    # data = np.genfromtxt('output/COMBINED_combined_snpy_cut_params.txt', dtype=str, delimiter=', ', skip_header=1)
    # dr3 = np.genfromtxt('txts/DR3_fits.dat', dtype=str, skip_header=1)
    # param_hist_compare(data[:, 15].astype(float), dr3[:, 25].astype(float),
    #                    0.1, title='91bg vs. DR3 -- EBVhost')
    # param_hist_compare(data[:,11].astype(float), dr3[:, 1].astype(float),
    #                    0.05, title='91bg vs. DR3 -- Stretch')

    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')
