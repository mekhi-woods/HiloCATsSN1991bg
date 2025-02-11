import os
import json
import glob
import utils  # Import of utils.py
import astropy
import requests
import numpy as np
import time as systime
import snpy as snoopyfit
import sncosmo as salt3fit
import matplotlib.pyplot as plt
from astropy.table import Table
from collections import OrderedDict
from astropy.time import Time as astrotime

class sneObj:
    # Construction Functions ----------------------------------------------------------------------------------------- #
    def __init__(self, source: str, algo: str, path: str):
        source = source.lower()
        if source == 'empty' or len(path) == 0:
            print('[+++] Creating empty SNe object...')
            pass
        elif source == 'class' and os.path.exists(path):
            print('[+++] Creating SNe object from class file...')
            self.load_class(path)
        elif source in ['atlas-91bg','atlas-norm', 'csp-91bg', 'csp-norm', 'ztf-91bg', 'ztf-norm'] and os.path.exists(path):
            print(f"[+++] Creating '{source}' SNe object using '{path}'...")

            # Extra details
            self.origin = source
            self.originalname = path.split('/')[-1].split('.')[0]
            self.path = f"classes/{source}/{source}_{self.originalname}_{algo}_class.txt"

            # Load and clean data
            var_tbl = self.make_objTbl(source, path)
            var_tbl = self.clean_objTbl(var_tbl)
            self.get_details(np.average(var_tbl['ra']), np.average(var_tbl['dec']))

            # Set arrays
            self.mjd = np.array(var_tbl['mjd'])
            self.mag = np.array(var_tbl['mag'])
            self.dmag = np.array(var_tbl['dmag'])
            self.flux = np.array(var_tbl['flux'])
            self.dflux = np.array(var_tbl['dflux'])
            self.filter = np.array(var_tbl['filter'])
            self.zp = np.array(var_tbl['zp'])
            self.params = {} # Initalized for later
            self.covariance = [] # Initalized for later

            # Save class file
            self.save_class()
        return
    def __str__(self):
        return (f"{self.objname}, {self.origin} @ {self.path}\n"
                f"{str(self.coords)} | ({self.z}, {str(round(self.z_cmb, 2))})")
    def make_objTbl(self, source: str, path: str) -> Table:
        # Load data and make intial table
        if source == 'atlas-91bg' or source == 'atlas-norm':
            # Load data
            with open(path, 'r') as f: hdr = f.readline()[1:-1].split(',')
            data = np.genfromtxt(path, dtype='str', delimiter=',', skip_header=1)

            # Make table
            var_table = Table()
            for h in hdr:
                try: var_table[h] = data[:, hdr.index(h)].astype(float)
                except ValueError: var_table[h] = data[:, hdr.index(h)]

            # Add parity with CSP & ZTF
            var_table.remove_columns(['err', 'chi/N', 'x', 'y', 'maj', 'min', 'phi', 'apfit', 'mag5sig', 'Sky', 'Obs'])
            var_table['zp'] = np.full(len(var_table), np.nan)
            for h_old, h_new in zip(['JD', 'm', 'dm', 'uJy', 'duJy', 'F', 'RA', 'Dec'],
                                    ['mjd', 'mag', 'dmag', 'flux', 'dflux', 'filter', 'ra', 'dec']):
                var_table[h_old].name = h_new
        elif source == 'csp-91bg' or source == 'csp-norm':
            var_table = Table(names=['filter', 'zp', 'z', 'ra', 'dec', 'mjd', 'mag', 'dmag', 'flux', 'dflux'],
                              dtype=[str, float, float, float, float, float, float, float, float, float])
            with open(path, 'r') as f:
                csp_objname, csp_z, csp_ra, csp_dec = f.readline()[2:-1].split(' ')
                for l in f.readlines():
                    l = l.split(' ')
                    # Filter line
                    if len(l) == 2:
                        csp_filter = str(l[1][:-1])
                        csp_zp = float(utils.get_constants()['csp_zpts_'+csp_filter])
                    else:
                        csp_mjd, csp_mag, csp_dmag = float(l[-3])+53000, float(l[-2]), float(l[-1])
                        csp_flux = 10 ** ((csp_mag - csp_zp) / -2.5)
                        csp_dflux = np.abs(csp_flux) * np.log(10) * ((1 / 2.5) * csp_dmag)

                        var_table.add_row([csp_filter, csp_zp, csp_z, csp_ra, csp_dec, csp_mjd,
                                           csp_mag, csp_dmag, csp_flux, csp_dflux])
            self.z = var_table['z'][0]  # Special case for CSP since CSP-I didn't report all SNe redshifts to TNS
        elif source == 'ztf-91bg' or source == 'ztf-norm':
            # Load hdr & data
            with open(path, 'r') as f:
                for i in range(3): f.readline()
                ztf_ra = float(f.readline().split(' ')[-2])
                ztf_dec = float(f.readline().split(' ')[-2])
            data = np.genfromtxt(path, delimiter=' ', dtype=str, skip_header=54)
            hdr, data = list(data[0]), data[1:]
            for i in range(len(hdr)): hdr[i] = hdr[i][:-1]

            # Make table
            var_table = Table()
            for h in hdr:
                try:
                    var_table[h] = data[:, hdr.index(h)].astype(float)
                except ValueError:
                    var_table[h] = data[:, hdr.index(h)]

            # Add RA & DEC
            var_table['ra'] = np.full(len(var_table), ztf_ra)
            var_table['dec'] = np.full(len(var_table), ztf_dec)

            # Fix time, JD to MJD
            var_table['jd'] = np.array(var_table['jd']).astype(float) - 2400000.5
            var_table['forcediffimflux'][var_table['forcediffimflux'] == 'null'] = 'nan'
            var_table['forcediffimfluxunc'][var_table['forcediffimfluxunc'] == 'null'] = 'nan'
            var_table['mag'] = (-2.5 * np.log10(np.array(var_table['forcediffimflux']).astype(float))) + var_table['zpdiff']
            var_table['dmag'] = np.abs(-1.08573620476 * (np.array(var_table['forcediffimfluxunc']).astype(float)
                                                         / np.array(var_table['forcediffimflux']).astype(float)))

            # Add parity with CSP & ATLAS
            var_table.remove_columns(['index', 'field', 'ccdid', 'qid', 'pid', 'infobitssci', 'sciinpseeing',
                                      'scibckgnd', 'scisigpix', 'zpmaginpsci', 'zpmaginpsciunc', 'zpmaginpscirms',
                                      'clrcoeff', 'clrcoeffunc', 'ncalmatches', 'exptime', 'adpctdif1', 'adpctdif2',
                                      'diffmaglim', 'programid', 'rfid', 'forcediffimsnr', 'forcediffimchisq',
                                      'forcediffimfluxap', 'forcediffimfluxuncap', 'forcediffimsnrap', 'aperturecorr',
                                      'dnearestrefsrc', 'nearestrefmag', 'nearestrefmagunc', 'nearestrefchi',
                                      'nearestrefsharp', 'refjdstart', 'refjdend', 'procstatu'])
            for h_old, h_new in zip(['zpdiff', 'jd', 'forcediffimflux', 'forcediffimfluxunc'],
                                    ['zp', 'mjd', 'flux', 'dflux']):
                var_table[h_old].name = h_new
        else:
            raise ValueError(f'[!!!] Unknown source, {source}! Must be [atlas-91bg/atlas-norm/csp/ztf]...')
        return var_table
    def clean_objTbl(self, tbl: Table) -> Table:
        # Remove '>' from ATLAS mags
        for i, n_mag in enumerate(tbl['mag']):
            if str(n_mag)[0] == '>': tbl['mag'][i] = float(n_mag[1:])

        # Remove nulls
        for col in ['mag', 'dmag', 'flux', 'dflux']:
            tbl = tbl[~np.isnan(np.array(tbl[col]).astype(float))]

        # Make cuts on magnitude errors
        mag_err_lim = float(utils.get_constants()['mag_err_lim'])
        if mag_err_lim > 0:
            tbl = tbl[np.array(tbl['dmag']).astype(float) < mag_err_lim]

        return tbl
    def get_details(self, ra: float, dec: float, check_key: bool = True):
        success, r = False, 2

        # Check TNS key first
        TNS_key_tb = utils.check_tnskey(ra, dec)
        if TNS_key_tb is not None and len(TNS_key_tb) == 1:
            print(f"[+++] Info found in TNS key! Pulling...")
            self.z_cmb = np.nan
            self.objname = TNS_key_tb['objname'][0]
            self.z = TNS_key_tb['z'][0]
            self.coords = [TNS_key_tb['ra'][0], TNS_key_tb['dec'][0]]
            self.discdate = TNS_key_tb['discoverydate'][0]
            return

        print(f"[+++] Querying TNS @ ({ra}, {dec}), with radius={r}...")
        while not success:
            try:
                details = get_TNSDetails(ra, dec, radius=r)

                # Special case for CSP since CSP-I didn't report all SNe redshifts to TNS
                if details['redshift'] is None and self.origin.split('-')[0] == 'csp':
                    print(f'[~~~] Speical CSP case: Keeping CSP heliocentric redshift {self.z}')
                    self.z = self.z
                else:
                    self.z = details['redshift']

                self.objname = details['objname']
                self.z_cmb = np.nan
                self.coords = [details['radeg'], details['decdeg']]
                self.discdate = float(astrotime(details['discoverydate'], format='iso').jd) - 2400000.5  # JD to MJD
                success = True
            except:
                print(f"[~~~] Warning: TNS found no targets @ ({ra}, {dec}), with radius={r}... Raising to r={r+1}")
                r += 1
            if r == 10:
                # raise RuntimeError('[!!!!!] TNS could not find transient!')
                print("[!!!!!] TNS could not find transient!")
                self.z_cmb = np.nan
                self.objname = self.originalname
                self.z = -99999
                self.coords = [-99999, -99999]
                self.discdate = -99999
                return None
        return
    def save_class(self):
        print(f"[+++] Saving {self.objname} class to {self.path}...")
        with open(self.path, 'w') as f:
            f.write(f"{self.origin},{self.objname},{self.originalname},{self.coords[0]},{self.coords[1]},{self.z},{self.z_cmb},{self.discdate}\n")
            f.write('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
            for p in self.params:
                f.write(f"{p},{self.params[p]['value']},{self.params[p]['err']}\n")
            f.write('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
            for line in self.covariance:
                p_line = str(line[0])
                for i in range(1, len(line)):
                    p_line += ',' + str(line[i])
                f.write(p_line + '\n')
            f.write('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
            for col, name in zip([self.zp, self.filter, self.mjd, self.flux, self.dflux, self.mag, self.dmag],
                                 ['zp', 'filter', 'mjd', 'flux', 'dflux', 'mag', 'dmag']):
                f.write(f"{name}")
                for n_col in col:
                    f.write(f",{n_col}")
                f.write('\n')
    def load_class(self, path: str):
        print(f"[+++] Opening class from {path}...")
        with open(path, 'r') as f:
            # Read header
            details = f.readline().split(',')
            self.path = path
            self.origin = details[0]
            self.objname = details[1]
            self.originalname = details[2]
            self.coords = [float(details[3]), float(details[4])]
            self.z = None if details[5] == "None" else float(details[5])
            self.z_cmb = float(details[6])
            self.discdate = float(details[7][:-1])
            f.readline()  # Skip break line

            # Read params
            self.params = {}
            line = f.readline()
            while '+++' not in line:
                line = line.split(',')
                self.params.update({line[0]: {'value': float(line[1]), 'err': float(line[2])}})
                line = f.readline()

            # Read covariances
            if 'salt' in path:
                self.covariance = np.array([])
                line = f.readline()
                while '+++' not in line:
                    line = line.split(',')
                    line[-1] = line[-1][:-1]
                    self.covariance = np.append(self.covariance, np.array(line).astype(float))
                    line = f.readline()
                if len(self.covariance) > 0:
                    self.covariance = self.covariance.reshape(4, 4)
            else:
                f.readline()  # Skip break line

            # Read arrays
            try:
                self.zp = np.array(f.readline().split(',')[1:])
                self.zp[-1] = self.zp[-1][:-1]
                self.filter = np.array(f.readline().split(',')[1:])
                self.filter[-1] = self.filter[-1][:-1]
                self.mjd = np.array(f.readline().split(',')[1:])
                self.mjd[-1] = self.mjd[-1][:-1]
                self.flux = np.array(f.readline().split(',')[1:])
                self.flux[-1] = self.flux[-1][:-1]
                self.dflux = np.array(f.readline().split(',')[1:])
                self.dflux[-1] = self.dflux[-1][:-1]
                self.mag = np.array(f.readline().split(',')[1:])
                self.mag[-1] = self.mag[-1][:-1]
                self.dmag = np.array(f.readline().split(',')[1:])
                self.dmag[-1] = self.dmag[-1][:-1]
            except IndexError:
                print(f"[~~~] While loading {self.objname}, empty arrays were found!")
                self.zp = np.array([])
                self.filter = np.array([])
                self.mjd = np.array([])
                self.flux = np.array([])
                self.dflux = np.array([])
                self.mag = np.array([])
                self.dmag = np.array([])
            except Exception as e:
                print(e)
        return

    # Display Functions ---------------------------------------------------------------------------------------------- #
    def simple_plot(self):
        for f in np.unique(self.filter):
            if f in ['J', 'H', 'Y', 'Ydw', 'u']: continue
            plt.scatter(self.mjd[self.filter == f].astype(float), self.flux[self.filter == f].astype(float), label=f)
        plt.axvline(self.discdate)
        plt.legend()
        plt.show()
        return
    def plot(self, y_type: str = 'mag', save_loc: str = ''):
        print(f"[+++] Plotting LC of '{self.objname}' ({y_type})...")
        filter_dict = {'u': 'teal', 'g': 'green', 'r': 'red', 'i': 'indigo', 'B': 'blue',
                       'V0': 'violet', 'V1': 'purple', 'V': 'red', 'Y': 'goldenrod', 'Hdw': 'tomato', 'H': 'salmon',
                       'J': 'aquamarine', 'Jrc2': 'cadetblue', 'Jdw': 'turquoise', 'Ydw': 'olive',
                       'c': 'cyan', 'o': 'orange', 'ZTF_g': 'green', 'ZTF_r': 'red', 'ZTF_i': 'indigo'}

        # Select y-axis element
        y_axis = self.mag if y_type == 'mag' else self.flux
        y_axis_err = self.dmag if y_type == 'mag' else self.dflux

        # Plot
        fig, axs = plt.subplots(1, 1, figsize=(12, 6), constrained_layout=True)
        for f in np.unique(self.filter):
            axs.errorbar(self.mjd[self.filter == f].astype(float), y_axis[self.filter == f].astype(float),
                         yerr=y_axis_err[self.filter == f].astype(float),
                         color = filter_dict[f], fmt='o', ms=3)

        # Lines
        axs.axvline(self.discdate, color='black', linestyle='--', label='Discovery Date')
        try: axs.axvline(self.params['Tmax']['value'], color='maroon', linestyle='--', label='Peak Brightness')
        except: pass

        # Formatting
        if y_type == 'mag': axs.invert_yaxis()
        axs.legend()
        axs.set_xlabel('MJD', size=16)
        axs.set_ylabel('Magnitude', size=16) if y_type == 'mag' else axs.set_ylabel('Flux (uJy)', size=16)
        plt.suptitle(f"Lightcurve of '{self.objname}' ({y_type})", size=16)
        if len(save_loc) != 0:
            print('[+++] Saving to '+save_loc)
            plt.savefig(save_loc)
        plt.show()
        return

    # Fitting Functions ---------------------------------------------------------------------------------------------- #
    def write_snpy_ascii(self, save_loc: str):
        filter_dict = {'o': 'ATri', 'c': 'ATgr', 't': 'ATri2',
                       'ZTF_g': 'g', 'ZTF_r': 'r', 'ZTF_i': 'i',
                       'B': 'B', 'H': 'H', 'J': 'J', 'Jrc2': 'Jrc2', 'V': 'V', 'V0': 'V0', 'Y': 'Y', 'Ydw': 'Ydw',
                       'g': 'g', 'i': 'i', 'r': 'r', 'u': 'u'}
        with open(save_loc, 'w') as f:
            # Line 1 -- Objname, Helio-Z, RA, Dec (Ex. SN1981D 0.005871 50.65992 -37.23272)
            f.write(f"{self.objname} {self.z} {self.coords[0]} {self.coords[1]}\n")
            for f_w in np.unique(self.filter):
                f_indexes = np.where(self.filter == f_w)[0]
                f.write(f"filter {filter_dict[f_w]}\n")
                for i in f_indexes:
                    # filter photometry block -- Date (JD/MJD), mag, err (i.e. 674.8593 12.94 0.11)
                    f.write(f"{self.mjd[i]} {self.mag[i]} {self.dmag[i]}\n")
        print(f'[+++] Saved file to... {save_loc}')
        return
    def snpy_fit(self):
        ascii_path = 'fitting/snpy-asciis/' + self.objname + '_ascii.txt'
        model_path = 'fitting/snpy-models/' + self.objname + '_model.txt'
        plot_path = 'fitting/snpy-plots/' + self.objname + '_lc.png'
        param_names = ['mu', 'st', 'Tmax', 'EBVhost']
        snpy_param_names = ['DM', 'st', 'Tmax', 'EBVhost']
        show_plots = True

        # Make ascii file for SNooPy to read
        self.write_snpy_ascii(ascii_path)

        # Load Data
        try:
            n_s = snoopyfit.get_sn(ascii_path)
        except Exception as error:
            self.params.update({'mu': {'value': 0.00, 'err': 0.00}})
            print('[!!!] Failed to load ASCII file -- ', error)
            return
        n_s.choose_model('EBV_model2', stype='st')
        n_s.set_restbands()  # Auto pick appropriate rest-bands

        # Remove empty filters -- fix for 'ValueError: attempt to get argmin of an empty sequence'
        for class_filter in list(n_s.data.keys()):
            if len(n_s.data[class_filter].magnitude) == 0:
                del n_s.data[class_filter]
            elif self.origin.split('-')[0] == 'csp' and class_filter in ['u', 'Y', 'J', 'H', 'Jrc2', 'Ydw']:
                print('[~~~] Speical CSP case: Removing ' + class_filter + '...')
                del n_s.data[class_filter]

        # Fit with SNooPy -- gives 5 tries before failing
        for i in range(5):
            print(f"[{i+1}/5] Attempting to fit '{self.objname}' with {list(n_s.data.keys())} =======================")
            try:
                ##### I care recall why we did this in the first place ####
                # if self.origin.split('-')[0] == 'csp':
                #     initial_filters = []
                #     for fil in ['B', 'V', 'g']:
                #         if fil in list(n_s.data.keys()):
                #             initial_filters.append(fil)
                #     print(f'[~~~] Speical CSP case: Fitting as {initial_filters} -> remaining...')
                #
                #     n_s.fit(initial_filters, dokcorr=True, k_stretch=False, reset_kcorrs=True,
                #             **{'mangle': 1, 'calibration': 0})
                #     n_s.fit(bands=None, dokcorr=True, k_stretch=False, reset_kcorrs=True,
                #             **{'mangle': 1, 'calibration': 0})
                # else:
                n_s.fit(bands=None, dokcorr=True, k_stretch=False, reset_kcorrs=True, **{'mangle': 1, 'calibration': 0})
                n_s.save(model_path)

                # Save parameters
                for j in range(len(param_names)):
                    self.params.update({param_names[j]: {'value': n_s.parameters[snpy_param_names[j]],
                                                         'err': n_s.errors[snpy_param_names[j]]}})
                self.params.update({'chisquare': {'value': n_s.model.chisquare,
                                                  'err': n_s.model.rchisquare}})
                n_s.plot(outfile=plot_path)
                if show_plots:
                    plt.show()
                    systime.sleep(2)
                plt.close()
                print(f'[+++] Successfully fit {self.objname}!')
                break
            except Exception as error:
                if 'All weights for filter' and 'are zero.' in str(error):
                    print('[!!!] Weights for filter', str(error).split(' ')[4], 'are zero. Removing...')
                    del n_s.data[str(error).split(' ')[4]]
                elif str(error) == 'Error:  to solve for EBVhost, you need to fit more than one filter':
                    print('[!!!] To few filters to fit!')
                    print(f'[---] Could not fit {self.objname}!')
                    self.params.update({'mu': {'value': -124.0, 'err': -124.0}})
                    break
                else:
                    self.params.update({'mu': {'value': -404.0, 'err': -404.0}})
                    print(error)
                    print(f'[---] Could not fit {self.objname}!')
                    break
        self.save_class()
        return
    def salt_fit(self):
        show_plot = True
        plot_path = f"fitting/salt-plots/{self.objname}_lc.png"

        CONSTANTS = utils.get_constants()
        alpha, beta = float(CONSTANTS['salt_alpha']), float(CONSTANTS['salt_beta'])
        mB_const, M0 = float(CONSTANTS['salt_mB_const']), float(CONSTANTS['salt_absolute_mag'])
        try:
            # Make sure more than one filter
            if len(np.unique(self.filter)) < 1:
                self.params.update({'mu': {'value': -124.0, 'err': -124.0}})
                self.save_class()
                print(f"[!!!] Too few filters! Can not fit {self.objname}! {np.unique(self.filter)}")
                return

            # Fix filters
            filter_dict = {'u': 'cspu', 'g': 'cspg', 'r': 'cspr', 'i': 'cspi', 'B': 'cspB',
                           'V0': 'cspv3014', 'V1': 'cspv3009', 'V': 'cspv9844', 'Y': 'cspys',
                           'J': 'cspjs', 'Jrc2': 'cspjd', 'Jdw': 'cspjd', 'Ydw': 'cspyd', 'Hdw': 'csphd', 'H': 'csphs',
                           'c': 'atlasc', 'o': 'atlaso', 'ZTF_g': 'ztfg', 'ZTF_r': 'ztfr', 'ZTF_i': 'ztfi'}
            salt_time, salt_filters, salt_flux = np.array([]), np.array([]), np.array([])
            salt_dflux, salt_zp = np.array([]), np.array([])

            for i in range(len(self.filter)):
                if self.origin.split('-')[0] == 'csp' and self.filter[i] in ['u', 'Y', 'J', 'H', 'Jrc2', 'Ydw']:
                    continue
                salt_time = np.append(salt_time, self.mjd[i])
                salt_filters = np.append(salt_filters, filter_dict[self.filter[i]])
                salt_flux = np.append(salt_flux, self.flux[i])
                salt_dflux = np.append(salt_dflux, self.dflux[i])
                salt_zp = np.append(salt_zp, self.zp[i])
            print('[~~~] Speical CSP case:', np.unique(self.filter), '->', np.unique(salt_filters))

            data = Table([salt_time, salt_filters, salt_flux, salt_dflux, salt_zp, np.full(len(salt_time), 'ab')],
                         names=('time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'))

            # Create Model
            model = salt3fit.Model(source='salt3')

            # Fit data to model
            model.set(z=self.z, t0=self.discdate)  # set the model's redshift.
            result, fitted_model = salt3fit.fit_lc(data, model, ['t0', 'x0', 'x1', 'c'], bounds={'x1': (-5, 5)})

            # Save Parameters
            param_names = ['t0', 'x0', 'x1', 'c']
            for i in range(len(param_names)):
                self.params.update({param_names[i]: {'value': result.parameters[i+1],
                                                     'err': result.errors[param_names[i]]}})

            # Save Covariance
            self.covariance = result['covariance']

            # Calculate
            pho_mB = -2.5 * np.log10(self.params['x0']['value']) + mB_const
            pho_mB_err = np.abs(-2.5 * (self.params['x0']['err'] / (self.params['x0']['value'] * np.log(10))))

            mu = pho_mB + (alpha * self.params['x1']['value']) - (beta * self.params['c']['value']) - M0
            mu_err = np.sqrt(pho_mB_err ** 2 + (np.abs(alpha) * self.params['x1']['err']) ** 2 + (np.abs(beta) * self.params['c']['err']) ** 2)

            self.params.update({'mu': {'value': mu, 'err': mu_err}})

            # Plot data with fit
            salt3fit.plot_lc(data, model=fitted_model, errors=result.errors)
            plt.savefig(plot_path)
            if show_plot:
                plt.show()
                systime.sleep(2)
            plt.close()

            self.save_class()
            print(f'[+++] Successfully fit {self.objname}!')
        except Exception as error:
            if 'result is NaN for' in str(error):
                print(f"[!!!] SALT3 couldn't fit with current parameter selection! Returning nan for {self.objname}...")
                self.params.update({'mu': {'value': -107.0, 'err': -107.0}})
            else:
                print(error)
                print(f'[---] Could not fit {self.objname}!')
                self.params.update({'mu': {'value': -404.0, 'err': -404.0}})
            self.save_class()
            return
        return
def get_TNSDetails(ra: str, dec: str, radius: str = '2'):
    """
    :param ra: Right Ascension in degrees
    :param dec: Declination in degrees
    :param radius: Radius to search out for in arcseconds
    :return: dict of TNS details
    """
    APIKEY = utils.get_apikeys()
    tns_bot_id, tns_bot_name, tns_bot_api_key = APIKEY['tns_bot_id'], APIKEY['tns_bot_name'], APIKEY['tns_bot_api_key']
    tns_marker = (
        f'tns_marker{{"tns_id": "{int(tns_bot_id)}",'
        f'"type": "bot", "name": "{tns_bot_name}"}}'
    )
    headers = {"User-Agent": tns_marker}
    search_obj = [
        ("ra", str(ra)),
        ("dec", str(dec)),
        ("radius", str(radius)),
        ("units", "arcsec"),
        ("objname", ""),
        ("objname_exact_match", 0),
        ("internal_name", ""),
        ("internal_name_exact_match ", 0),
        ("objid", ""),
        ("public_timestamp", ""),
    ]
    search_data = {"api_key": tns_bot_api_key, "data": json.dumps(OrderedDict(search_obj))}
    response = requests.post("https://www.wis-tns.org/api/get/search", headers=headers, data=search_data)
    response = json.loads(response.text)
    transients =  response["data"]
    get_obj = [
        ("objname", transients[0]["objname"]),
        ("objid", transients[0]["objid"]),
        ("photometry", "0"),
        ("spectra", "0"),
    ]
    get_data = {"api_key": tns_bot_api_key, "data": json.dumps(OrderedDict(get_obj))}
    response = requests.post("https://www.wis-tns.org/api/get/object", headers=headers, data=get_data)
    response = json.loads(response.text)
    details = response["data"]
    return details
def fit_subprocess(dataset: str, path: str, algo: str, rewrite: bool = False):
    """
    :param dataset: Data set to pull light curves from
    :param path: Path to pull lightcurve
    :param algo: Algorithm to use for fitting
    :param rewrite: Even if data is already procced, it will act as if its the not
    :return: sneObj object
    """
    # Check if preveiously fit
    class_save_loc = f"classes/{dataset}/{dataset}_{path.split('/')[-1].split('.txt')[0]}_{algo}_class.txt"
    print(f"[+++] Checking for saved class at {class_save_loc}...")
    if not rewrite and os.path.exists(class_save_loc):
        print(f"[+++] Class found! Loading from {class_save_loc}...")
        sn = sneObj('class', algo, class_save_loc)
    else:
        print(f"[---] No class found! Setting up class light curve at {path}...")
        sn = sneObj(dataset, algo, path)

    # Fit SN class
    if len(sn.params) == 0 or rewrite:
        print(f"[+++] Fitting '{sn.objname}' with '{dataset}' data & the '{algo}' algorithm...")
        if algo.lower() == 'snpy':
            sn.snpy_fit()
        elif algo.lower() == 'salt':
            sn.salt_fit()
    else:
        print(f"[+++] '{sn.objname}' already fit! Loading...")

    update_TNS = True
    if update_TNS:
        if utils.check_tnskey(sn.coords[0], sn.coords[1]) is None:
            print(f"[+++] Updating TNS with '{sn.objname}' data...")
            utils.append_tnskey(sn.coords[0], sn.coords[1], sn.objname, sn.z, sn.discdate)
    return sn
def fit(data_loc: str, algo: str, rewrite: bool = False) -> list[sneObj]:
    """
    :param data_loc: Location of data; if single path -> indivisual mode, if directory -> batch mode
    :param algo: Algorithm to fit; either SNooPy or SALT3
    :param rewrite: Even if data is already procced, it will act as if its the not
    :return: sneObj object or list of sneObj objects
    """
    paths = glob.glob(data_loc)
    dataset = paths[0].split('/')[-2].lower()

    # Verify proper dataset & proper algorithm
    valid_datasets = ['csp-91bg', 'csp-norm', 'atlas-91bg','atlas-norm', 'ztf-91bg', 'ztf-norm']
    valid_algorithms = ['snpy', 'salt']
    if dataset not in valid_datasets: raise ValueError(f"[!!!] Dataset, '{dataset}', not recognized! {valid_datasets}")
    elif algo not in valid_algorithms: raise ValueError(f"[!!!] Algorithm, '{algo}', not recognized! {valid_algorithms}")

    # Select fitting mode
    ## Indivisual mode
    if len(paths) == 1:
        print(f"[+++] Fitting file '{data_loc}'...")  # Indivisual fit
        sn = fit_subprocess(dataset, paths[0], algo, rewrite)
        return [sn]
    ## Batch mode
    elif len(paths) > 1:
        print(f"[+++] Fitting data in '{data_loc}'...")  # Batch fit
        success_counter, sne, fail_reasons  = 0, [], []
        for i, path in enumerate(paths):
            print(f'[{i + 1} / {len(paths)}] ================================================================')
            sn = fit_subprocess(dataset, path, algo, rewrite)

            # Check if success
            if sn.discdate == -99999:  # TNS Faliure
                fail_reasons.append(f"{sn.objname}: TNS Faliure! Needs manual TNS")
            elif sn.params['mu']['value'] == -124.0:
                fail_reasons.append(f"{sn.objname}: Not enough filters to fit!")
            elif sn.params['mu']['value'] == -107.0:
                fail_reasons.append(f"{sn.objname}: SALT3 couldn't fit with current parameter selection!")
            elif sn.params['mu']['value'] < 0:
                fail_reasons.append(f"{sn.objname}: Unknown error!")
            else:
                success_counter += 1
                sne.append(sn)
        print(f"[!!!!!] Successfully fit {success_counter} / {len(paths)}!")

        # Show error
        show_errors = True
        if show_errors:
            for n in fail_reasons: print(n)  # Show failed fit reasons

        # Sorted error table for README.md
        readme_errors = True
        if readme_errors:
            print(f"{data_loc.split('/')[1].split('-')[0]}+{algo.upper()}+{data_loc.split('/')[1].split('-')[1]}")
            low_filters, tns_needed, unknown = [], [], []
            for i, n in enumerate(fail_reasons):
                if n.split(': ')[-1] == "Not enough filters to fit!": low_filters.append(n.split(': ')[0])
                elif n.split(': ')[-1] == "TNS Faliure! Needs manual TNS": tns_needed.append(n.split(': ')[0])
                else: unknown.append(n.split(': ')[0])
            print("Not enough filters to fit!")
            for i, n in enumerate(low_filters):
                if (i+1) % 5 == 0:
                    print("SN" + n, end=",<br/> ")
                else:
                    print("SN"+n, end=", ")
            print("\nTNS Faliure! Needs manual TNS!")
            for i, n in enumerate(tns_needed):
                if (i+1) % 5 == 0:
                    print("SN" + n, end=",<br/> ")
                else:
                    print("SN"+n, end=", ")
            print('\n')
            print("\nUnknown error!")
            for i, n in enumerate(unknown):
                if (i+1) % 5 == 0:
                    print("SN" + n, end=",<br/> ")
                else:
                    print("SN"+n, end=", ")
            print('\n')
        return sne
    else:
        print('[!!!] Invalid file/data path!')
    return


if __name__ == '__main__':
    start = systime.time()  # Runtime tracker
    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')
