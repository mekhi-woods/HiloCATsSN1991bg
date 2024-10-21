# M.D. Woods - 10/17/24
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
from scripts import get_vpec  # Loads heavy data so if not activly fitting I would comment this out
from scripts.salt3_param_fitter import optimize_alpha_beta

CONSTANTS = gen.get_constants()

class SN91bg:
    """
    A class to store file data for CSP/ATLAS/ZTF and fit data of both SNooPy & SALT3.
    """
    def __init__(self, data_set: str = 'EMPTY', path: str = '', dmag_max: float = 0.00, dflux_max: float = 0.00):
        """
        The initialization function for the class. Takes a file of data and preps it for fitting.
        :param data_set: str = 'EMPTY'; name of the data set
        :param path: str; location of the data file
        :param dmag_max: float = 0.00; maximum desired magnitude value
        :param dflux_max: float = 0.00; maximum desired flux value
        """
        # Check if file empty
        if len(path) != 0:
            with open(path, 'r') as f:
                if len(f.readlines()) < 3:
                    print('[!!!] File [' + path + '] empty!')
                    data_set = 'EMPTY'

        # If no data set given, try to guess
        if len(data_set) == 0:
            if 'CSP' in path.upper():
                data_set = 'CSP'
            elif 'ATLAS' in path.upper():
                data_set = 'ATLAS'
            elif 'ZTF' in path.upper():
                data_set = 'ZTF'
            else:
                raise ValueError('[!!!] Data set not recognized! [CSP/ATLAS/ZTF]')

        # Unpack data file based on data set origin
        if data_set.upper() == 'CSP':
            print('[+++] Creating class using CSP data...')
            print(path)

            # Set static variables
            self.origin = 'CSP'
            self.z_cmb = np.nan
            self.params = {}
            self.covariance = np.array([])

            # Header elements -- receive original, z, ra, dec
            with open(path, 'r') as f:
                hdr = f.readline().split(' ')
                self.originalname, self.z = hdr[0], float(hdr[1])
                self.coords = [float(hdr[2]), float(hdr[3][:-1])]

            # Query TNS for transient details -- receive objname & discovery date
            self.objname, z_void, self.discovery_date = gen.TNS_details(self.coords[0], self.coords[1])  # Void TNS redshift

            # Initialize empty arrays
            self.zp, self.filters, self.time, self.flux, self.dflux, self.mag, self.dmag = (
                np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]))

            # Set arrays
            n_filter = ''
            with open(path, 'r') as f:
                f.readline()  # Skips header
                for line in f.readlines():
                    data_line = line[:-1].split(' ')
                    if len(data_line) == 2:
                        n_filter = data_line[1]
                    elif len(data_line) >= 3:
                        if len(data_line) == 4:
                            data_line = data_line[1:]
                        n_time, n_mag, n_dmag = float(data_line[0]), float(data_line[1]), float(data_line[2])
                        n_zp = float(CONSTANTS['csp_zpts_' + n_filter])
                        n_time = n_time + 53000  # JD to MJD

                        n_flux = 10 ** ((n_mag - n_zp) / -2.5)
                        n_dflux = np.abs(n_flux) * np.log(10) * ((1 / 2.5) * n_dmag)

                        self.zp = np.append(self.zp, n_zp)
                        self.filters = np.append(self.filters, n_filter)
                        self.time = np.append(self.time, n_time)
                        self.flux = np.append(self.flux, n_flux)
                        self.dflux = np.append(self.dflux, n_dflux)
                        self.mag = np.append(self.mag, n_mag)
                        self.dmag = np.append(self.dmag, n_dmag)

            self.period = [np.min(self.time), np.max(self.time)]
        elif data_set.upper() == 'ATLAS':
            print('[+++] Creating class using ATLAS data...')
            print(path)

            # Set static variables
            self.origin = 'ATLAS'
            self.z_cmb = np.nan
            self.params = {}
            self.covariance = np.array([])

            # Get original name from file name
            self.originalname = path.split('/')[-1][:-4]

            # Load data
            data = np.genfromtxt(path, delimiter=',', dtype=str, skip_header=1)

            # Query TNS for transient details
            self.coords = [np.average(data[:, 1].astype(float)), np.average(data[:, 2].astype(float))]
            self.objname, self.z, self.discovery_date = gen.TNS_details(self.coords[0], self.coords[1])
            self.discovery_date = float(self.discovery_date)
            if self.z == 'None': self.z = np.nan
            else: self.z = float(self.z)

            # Set arrays
            self.zp = data[:, 7].astype(float)
            self.filters = data[:, 6]
            self.time = data[:, 8]
            self.flux = data[:, 16]
            self.dflux = data[:, 17]
            self.mag = data[:, 3]
            self.dmag = data[:, 4]

            # Remove '<' & '>' from magnitudes
            for n in range(len(self.mag)):
                self.mag[n] = self.mag[n].replace('>', '')
                self.mag[n] = self.mag[n].replace('<', '')
        elif data_set.upper() == 'ZTF':
            print('[+++] Creating class using ZTF data...')
            print(path)

            # Set static variables
            self.origin = 'ZTF'
            self.z_cmb = np.nan
            self.params = {}
            self.covariance = np.array([])

            # Load data
            data = np.genfromtxt(path, delimiter=None, dtype=str, skip_header=56)

            # Get original name from file name
            self.originalname = path.split('/')[-1].split('.')[0].split('_')[1]

            with (open(path, 'r') as f):
                # Open header
                hdr = f.readlines()
                self.coords = [float(hdr[3].split(' ')[-2]), float(hdr[4].split(' ')[-2])]
                self.objname, self.z, self.discovery_date = gen.TNS_details(self.coords[0], self.coords[1])
                self.z = np.nan if self.z == 'None' else float(self.z) # Remove?

                # Initialize arrays
                self.time, self.flux, self.dflux, self.zp, self.filters = (
                    data[:, 22], data[:, 24], data[:, 25], data[:, 20], data[:, 4])

                # Remove nulls from arrays
                self.time = self.time[(self.flux != 'null')].astype(float)
                self.zp = self.zp[(self.flux != 'null')].astype(float)
                self.filters = self.filters[(self.flux != 'null')]
                self.dflux = self.dflux[(self.flux != 'null')].astype(float)
                self.flux = self.flux[(self.flux != 'null')].astype(float)

                # Fix time, JD to MJD
                self.time = self.time - 2400000.5

                # Flux -> Mag
                self.mag = (-2.5 * np.log10(self.flux)) + self.zp  # Get magnitudes m = -2.5log(F) + zp
                self.dmag = np.abs(-1.08573620476 * (self.dflux / self.flux))

            # Adjusting around tmax
            ztf_spread = float(CONSTANTS['ztf_spread'])
            self.discovery_date = float(self.discovery_date)
            if ztf_spread != 0 and len(self.time) != 0:
                self.flux = self.flux[(self.time > self.discovery_date - ztf_spread) &
                                      (self.time < self.discovery_date + ztf_spread)]
                self.dflux = self.dflux[(self.time > self.discovery_date - ztf_spread) &
                                        (self.time < self.discovery_date + ztf_spread)]
                self.mag = self.mag[(self.time > self.discovery_date - ztf_spread) &
                                    (self.time < self.discovery_date + ztf_spread)]
                self.dmag = self.dmag[(self.time > self.discovery_date - ztf_spread) &
                                      (self.time < self.discovery_date + ztf_spread)]
                self.zp = self.zp[(self.time > self.discovery_date - ztf_spread) &
                                  (self.time < self.discovery_date + ztf_spread)]
                self.filters = self.filters[(self.time > self.discovery_date - ztf_spread) &
                                            (self.time < self.discovery_date + ztf_spread)]
                self.time = self.time[(self.time > self.discovery_date - ztf_spread) &
                                      (self.time < self.discovery_date + ztf_spread)]
        elif data_set.upper() == 'EMPTY' or len(path) == 0:
            print('[+++] Creating an empty class...')
            self.objname = ''
            self.originalname = ''
            self.coords = [np.nan, np.nan]
            self.z = np.nan
            self.z_cmb = np.nan
            self.origin = ''
            self.discovery_date = np.nan

            self.period = [np.nan, np.nan]
            self.params = {}
            self.covariance = np.array([])

            self.zp = np.array([])
            self.filters = np.array([])
            self.time = np.array([])
            self.flux = np.array([])
            self.dflux = np.array([])
            self.mag = np.array([])
            self.dmag = np.array([])
        else:
            raise ValueError("Data set '" +
                             data_set + "' not recognized [CSP/ATLAS/ZTF/EMPTY]")

        # Remove None
        indexes = (self.mag != 'None') & (self.flux != 'None')
        self.time = self.time[indexes].astype(float)
        self.zp = self.zp[indexes].astype(float)
        self.filters = self.filters[indexes]
        self.dflux = self.dflux[indexes].astype(float)
        self.flux = self.flux[indexes].astype(float)
        self.dmag = self.dmag[indexes].astype(float)
        self.mag = self.mag[indexes].astype(float)

        # Remove zero errors
        indexes = (self.dmag > 0.00) & (self.flux > 0.00)
        self.time = self.time[indexes]
        self.zp = self.zp[indexes]
        self.filters = self.filters[indexes]
        self.dflux = self.dflux[indexes]
        self.flux = self.flux[indexes]
        self.dmag = self.dmag[indexes]
        self.mag = self.mag[indexes]

        print('[+++] ' + self.objname + ' -- Cleaning data...')
        # Adjust Maximums
        if dmag_max == 'median': dmag_max = np.median(self.dmag)
        if dflux_max == 'median': dflux_max = np.median(self.dflux)
        if dmag_max == 'average': dmag_max = np.median(self.dmag)
        if dflux_max == 'average': dflux_max = np.average(self.dflux)

        # Cut arrays based on dmag_max/dflux_max
        if dmag_max != 0.00:
            self.flux = self.flux[self.dmag < dmag_max]
            self.dflux = self.dflux[self.dmag < dmag_max]
            self.filters = self.filters[self.dmag < dmag_max]
            self.zp = self.zp[self.dmag < dmag_max]
            self.time = self.time[self.dmag < dmag_max]
            self.mag = self.mag[self.dmag < dmag_max]
            self.dmag = self.dmag[self.dmag < dmag_max]
        if dflux_max != 0.00:
            self.mag = self.mag[self.dflux < dflux_max]
            self.dmag = self.dmag[self.dflux < dflux_max]
            self.filters = self.filters[self.dflux < dflux_max]
            self.zp = self.zp[self.dflux < dflux_max]
            self.time = self.time[self.dflux < dflux_max]
            self.flux = self.flux[self.dflux < dflux_max]
            self.dflux = self.dflux[self.dflux < dflux_max]

        # Final stage
        self.period = None
        if len(self.time) > 0:
            self.period = [np.min(self.time), np.max(self.time)]
            print('[+++] Class created successfully!')
            return
        else:
            print('[!!!] No valid points found in file!')
            return
    # Display Options ----------------------------------------------------------------------------------------------- #
    def __str__(self):
        """
        Displays data in single line.
        """
        return (self.objname + ' | ' + self.originalname + ' | ' + self.origin + ' | ' + str(self.coords) +
                ' | (z | z_cmb): (' + str(self.z) + ' | ' + str(round(self.z_cmb, 2)) + ')')
    def print_info(self):
        """
        Displays important data in readable format.
        """
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
        return
    def plot(self, y_type: str = 'mag', save_loc: str = '', zoom: int = 0, subplots: bool = False, date_lines: bool = True):
        """
        Plots the SNe lightcurve.
        :param y_type: str = 'mag'; either mag or flux, the unit to plot on the y-axis
        :param save_loc: str = ''; place to save the plot
        :param zoom: int = 0; number of MJD to zoom into around tmax
        :param subplots: bool = False; either plot each filter together or indivisual
        :param date_lines: bool = True; show tmax and discovery date as vertical line on plot
        """
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

    # Loading/Saving Data ------------------------------------------------------------------------------------------- #
    def save_class(self, save_loc: str):
        """
        Saves class to file
        :param save_loc: str; location to save class to
        """
        print('[+++] '+self.objname+' -- Saving class to '+save_loc+'classes/'+self.objname+'_class.txt')
        with open(save_loc+'classes/'+self.objname+'_class.txt', 'w') as f:
            f.write(self.origin + ' ' + self.objname + ' ' + self.originalname +
                    ' ' + str(self.coords[0]) + ' ' + str(self.coords[1]) +
                    ' ' + str(self.z) + ' ' + str(self.z_cmb) + ' ' + str(self.discovery_date) + '\n')

            f.write('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
            for p in self.params:
                f.write(p+', ' + str(self.params[p]['value']) + ', ' + str(self.params[p]['err']) + '\n')

            f.write('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
            for line in self.covariance:
                p_line = str(line[0])
                for i in range(1, len(line)):
                    p_line += ', ' + str(line[i])
                f.write(p_line+'\n')

            f.write('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
            cats = [self.zp, self.filters, self.time, self.flux, self.dflux, self.mag, self.dmag]
            cat_names = ['zp', 'filters', 'time', 'flux', 'dflux', 'mag', 'dmag']
            for i in range(len(cat_names)):
                f.write(cat_names[i])
                for j in range(len(self.zp)):
                    f.write(', ' + str(cats[i][j]))
                f.write('\n')
        return
    def load_from_file(self, file: str):
        """
        Loads class from file
        :param file: str; location to restore class from
        """
        print('[+++] Loading class from '+file)
        with open(file, 'r') as f:
            # Get header data
            hdr = f.readline().split(' ')
            self.origin, self.objname, self.originalname, self.coords, self.z, self.z_cmb, self.discovery_date = (
                hdr[0], hdr[1], hdr[2],
                (float(hdr[3]), float(hdr[4])), float(hdr[5]), float(hdr[6]), float(hdr[7][:-1]))

            f.readline() # Skip +++++++ line

            # Read in parameter data
            line = f.readline()
            while '+++' not in line:
                line = line.split(', ')
                self.params.update({line[0]: {'value': float(line[1]), 'err': float(line[2])}})
                line = f.readline()

            # Load covariance from file
            if 'salt' in file:
                self.covariance = np.array([])
                line = f.readline()
                while '+++' not in line:
                    line = line.split(', ')
                    line[-1] = line[-1][:-1]
                    self.covariance = np.append(self.covariance, np.array(line).astype(float))
                    line = f.readline()
                self.covariance = self.covariance.reshape(4, 4)

            # Load data arrays
            self.zp = np.array(f.readline().split(', ')[1:]).astype(float)
            self.filters = np.array(f.readline().split(', ')[1:]).astype(str)
            self.filters[-1] = self.filters[-1][:-1]  # Removes the /n from the end of the line
            self.time = np.array(f.readline().split(', ')[1:]).astype(float)
            self.flux = np.array(f.readline().split(', ')[1:]).astype(float)
            self.dflux = np.array(f.readline().split(', ')[1:]).astype(float)
            self.mag = np.array(f.readline().split(', ')[1:]).astype(float)
            self.dmag = np.array(f.readline().split(', ')[1:]).astype(float)

            self.period = (np.min(self.time), np.max(self.time))
        return

    # Fitting Functions --------------------------------------------------------------------------------------------- #
    def fit(self, algo: str = 'snpy', save_loc: str = ''):
        """
        Selects the algorithm to fit class data with
        :param algo: str; algorithm to use, either snpy or salt
        :param save_loc: str; location to save class to
        """
        if len(algo) == 0: return False
        if len(save_loc) == 0: save_loc = CONSTANTS[algo+'_'+self.origin.lower() + '_saved_loc']

        # Choose algorithm
        if algo == "snpy":
            self.write_snpy_ascii(save_loc=save_loc + 'ascii/')
            self.snpy_fit(save_loc=save_loc)
        elif algo == 'salt':
            self.salt_fit(save_loc=save_loc)

        # Check for successful fit
        if (self.params['mu']['value'] <= 0.00):
            print(f"[!!!] Invald distance modlus: {self.params['mu']['value']}")
            return False

        # Get Host Mass
        self.get_host_mass(use_key=True)

        # Check for vaild mass
        if (self.params['hostMass']['value'] <= 0.00):
            print(f"[!!!] Invald host mass: {self.params['hostMass']['value']}")
            return False

        # Save
        self.save_class(save_loc)
        return True
    def snpy_fit(self, save_loc: str, use_saved: bool = False, show_plot: bool = True, quiet: bool = False):
        """
        Fit class data using snpy algorithm
        :param save_loc: str; location to save class to
        :param use_saved: bool = False; toggled to use a saved class
        :param show_plot: bool = True; toggled to show plot
        :param quiet: bool = False; toggled to turn on quiet mode
        """
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

        # Fit with SNooPy -- gives 5 tries before failing
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
                self.params.update({'chisquare': {'value': n_s.model.chisquare,
                                                  'err': n_s.model.rchisquare}})

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
                    self.params.update({'mu': {'value': -999.0, 'err': -999.0}})
                    break
                else:
                    self.params.update({'mu': {'value': -1.0, 'err': -1.0}})
                    print(error)

        # Restore print statements
        sys.stdout = sys.__stdout__

        print('[+++] Successfully fit ' + self.objname + '!')
        return
    def salt_fit(self, save_loc: str, show_plot: bool = True, quiet: bool = False):
        """
        Fit class data using salt algorithm
        :param save_loc: str; location to save class after fitting
        :param show_plot: bool = True; toggle to show plot
        :param quiet: bool = False; toggle to quiet mode
        """
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
            model.set(z=self.z, t0=self.discovery_date)  # set the model's redshift.
            result, fitted_model = sncosmo.fit_lc(data, model, ['t0', 'x0', 'x1', 'c'], bounds={'x1': (-5, 5)})

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

            mu = pho_mB + (alpha * self.params['x1']['value']) - (beta * self.params['c']['value'])  - M0
            mu_err = np.sqrt(pho_mB_err ** 2 + (np.abs(alpha) * self.params['x1']['err']) ** 2 + (np.abs(beta) * self.params['c']['err']) ** 2)

            self.params.update({'mu': {'value': mu, 'err': mu_err}})

            # Plot data with fit
            if show_plot:
                sncosmo.plot_lc(data, model=fitted_model, errors=result.errors)
                print('Saving plots to', save_loc + '/plots/' + self.objname + '_salt3lc.png')
                plt.savefig(save_loc + '/plots/' + self.objname + '_salt3lc.png')
                plt.show()

            print('[+++] Successfully fit ' + self.objname + '!')
        except Exception as error:
            print(error)
            self.params.update({'mu': {'value': -1.0, 'err': -1.0}})

        # Restore print statements
        sys.stdout = sys.__stdout__

        return
    def write_snpy_ascii(self, save_loc: str):
        """
        Writes data to snpy ascii file.
        :param save_loc: str; location to save ascii file
        """
        print('[+++] '+self.objname+' -- Saving data to ASCII files for SNooPy...')
        filter_dict = {'o': 'ATri', 'c': 'ATgr', 't': 'ATri2',
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
    def get_host_mass(self, use_key: bool = False):
        """
        Gets host mass of SN using GHOST
        :param use_key: bool; toggle to use key to skip finding host mass
        """
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
                    self.params.update({'hostMass': {'value': -999, 'err': -999}})
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
def norm_fit(algo: str = 'snpy', save_loc: str = '', dmag_max: float = 0.00, dflux_max: float = 0.00) -> list[object]:
    """
    Top level function to easily fit normal SNe Ia with either snpy or salt
    :param algo: str = 'snpy'; algorithm to use for fitting
    :param save_loc: str = ''; location to save fits paramaters
    :param dmag_max: float = 0.00; maximum dmag to allow when loading data
    :param dflux_max: float = 0.00; maximum dflux to allow when loading data
    :return: SNe: list[object]; List of SN91bg classes from fitting call.
    """
    SNe, files = [], glob.glob('data/CSPdata/*.txt')
    for path in files:
        print('[', files.index(path) + 1, '/', len(files), ']')
        print('-----------------------------------------------------------------------------------------------')
        tempSN = SN91bg(path=path, data_set = 'CSP', dmag_max = dmag_max, dflux_max = dflux_max)
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
    if len(save_loc) != 0:
        save_params_to_file(save_loc, SNe)
    return SNe
def combined_fit(algo: str = 'snpy', dmag_max: float = 0.00, dflux_max: float = 0.00) -> list[object]:
    """
    Top level function to easily fit CSP/ATLAS/ZTF 91bg-like SNe Ia with either snpy or salt; accounts for overlap with
    ATLAS & ZTF
    :param algo: str = 'snpy'; algorithm to use for fitting
    :param dmag_max: float = 0.00; maximum dmag to allow when loading data
    :param dflux_max: float = 0.00; maximum dflux to allow when loading data
    :return: fit_combined_SNe: list[object]; List of SN91bg classes from fitting call.
    """
    sys.stdout = open(os.devnull, 'w')  # Lots of unnecessary output

    # Load and set up files
    csp_files, atlas_files, ztf_files = glob.glob('data/CSP/*.txt'), glob.glob('data/ATLAS/*.txt'), glob.glob('data/ZTF/*.txt')
    CSP_SNe, ATLAS_SNe, ZTF_SNe, atlas_names = {}, {}, {}, []
    for file in csp_files:
        tempSN = SN91bg(path=file, data_set = 'CSP', dmag_max = dmag_max, dflux_max = dflux_max)
        if len(tempSN.objname) != 0: CSP_SNe.update({tempSN.objname: tempSN})
    for file in atlas_files:
        tempSN = SN91bg(path=file, data_set = 'ATLAS', dmag_max = dmag_max, dflux_max = dflux_max)
        if len(tempSN.objname) != 0: ATLAS_SNe.update({tempSN.objname: tempSN})
    for file in ztf_files:
        tempSN = SN91bg(path=file, data_set = 'ZTF', dmag_max = dmag_max, dflux_max = dflux_max)
        if len(tempSN.objname) != 0: ZTF_SNe.update({tempSN.objname: tempSN})

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
        print('[', combined_SNe.index(n_SN) + 1, '/', len(combined_SNe), '] Fitting data for ' + n_SN.objname + ' [' + algo + '] [' + n_SN.origin + ']...')
        print('-----------------------------------------------------------------------------------------------')

        # Fit class
        successfulFit = n_SN.fit(algo)

        # Check if fit fail in any way
        if not successfulFit:
            print(f"[-----] An error occured while fitting! Discarding {n_SN.objname}...")
            continue
        else:
            print('[+++++] Success!')
            fit_combined_SNe.append(n_SN)

    print('=====================================================================================================\n',
          'Successfully fit [', len(fit_combined_SNe), '/', len(combined_SNe), ']!\n',
          '=====================================================================================================\n')
    return fit_combined_SNe
def batch_fit(data_set: str, algo: str = 'snpy', dmag_max: float = 0.00, dflux_max: float = 0.00) -> list[object]:
    """
    Top level function to easily batch fit 91bg-like SNe Ia with snpy/salt using CSP/ATLAS/ZTF
    :param data_set: str; choice of dataset [CSP/ATLAS/ZTF]
    :param algo: str = 'snpy'; algorithm to use for fitting
    :param dmag_max: float = 0.00; maximum dmag to allow when loading data
    :param dflux_max: float = 0.00; maximum dflux to allow when loading data
    :return: SNe: list[object]; List of SN91bg classes from fitting call.
    """
    SNe, files = [], glob.glob(CONSTANTS[data_set.lower()+'_data_loc'] + '*.txt')
    for path in files:
        print('[', files.index(path) + 1, '/', len(files), ']')
        print('-----------------------------------------------------------------------------------------------')
        tempSN = SN91bg(path=path, data_set = data_set, dmag_max = dmag_max, dflux_max = dflux_max)

        # Check if data file could make class
        if len(tempSN.objname) == 0: continue

        # Fit class
        successfulFit = tempSN.fit(algo)

        # Check if fit fail in any way
        if not successfulFit:
            print(f"[-----] An error occured while fitting! Discarding {tempSN.objname}...")
            continue
        else:
            print('[+++++] Success!')
            SNe.append(tempSN)

    print('Sucessfully fit [', len(SNe), '/', len(files), ']!')
    return SNe
def batch_load(data_set: str, algo: str = 'snpy') -> list[object]:
    """
    Loads class data based on dataset
    :param data_set: str; choice of dataset [CSP/ATLAS/ZTF]
    :param algo: str = 'snpy'; algorithm to use for fitting
    :return: SNe: list[object]; List of SN91bg classes from fitting call.
    """
    SNe = []
    for path in glob.glob('saved/' + algo + '/' + data_set.lower() + '/classes/*_class.txt'):
        sys.stdout = open(os.devnull, 'w')
        tempSN = SN91bg()
        sys.stdout = sys.__stdout__

        tempSN.load_from_file(path)
        SNe.append(tempSN)
    return SNe

# File Saving / Augmentation ---------------------------------------------------------------------------------------- #
def save_params_to_file(save_loc: str, SNe: list[object]):
    """
    Saves parameters of multiple SNe to single file
    :param save_loc: str; location to save paramter file
    :param SNe: list[object]; list of SNe objects
    """
    print('[+++] Saving params to '+save_loc+'...')

    # Check if you can save
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
def save_params_to_file_cov(save_loc: str, SNe: list[object]):
    """
    Saves parameters of multiple SNe to single file
    :param save_loc: str; location to save paramter file
    :param SNe: list[object]; list of SNe objects
    """
    print('[+++] Saving params to ' + save_loc + '...')

    # Check if you can save
    if len(SNe) == 0:
        print('[+++] No params to save!')
        return

    # Build header
    hdr = ['objname', 'ra', 'dec', 'z', 'z_cmb', 'MJDs', 'MJDe', 'origin']
    for cat in SNe[0].params.keys():
        hdr.append(cat)
        hdr.append(cat+'_err')
    if len(SNe[0].covariance) != 0:
        hdr = hdr + ['cov00', 'cov01', 'cov02', 'cov03',
                     'cov10', 'cov11', 'cov12', 'cov13',
                     'cov20', 'cov21', 'cov22', 'cov23',
                     'cov30', 'cov31', 'cov32', 'cov33']

    # Write data to file
    with open(save_loc, 'w') as f:
        # Write header
        f.write(str(hdr[0]))
        for i in range(1, len(hdr)):
            f.write(', '+hdr[i])
        f.write('\n')

        # Write data by SN
        for n_SN in SNe:
            # Write header variables
            line = (n_SN.objname + ', ' + str(n_SN.coords[0]) + ', ' + str(n_SN.coords[1]) + ', ' +
                    str(n_SN.z) + ', ' + str(n_SN.z_cmb) + ', ' +
                    str(n_SN.period[0]) + ', ' + str(n_SN.period[1]) + ', ' + str(n_SN.origin))

            # Write paramaters
            for param in n_SN.params:
                line += ', ' + str(n_SN.params[param]['value']) + ', ' + str(n_SN.params[param]['err'])

            # Write covariance
            if len(SNe[0].covariance) != 0:
                for arr_line in n_SN.covariance:
                    for cov in arr_line:
                        line += ', ' + str(cov)

            # Commit to file line
            f.write(line + '\n')
    return
def sample_cutter(path: str, algo: str ='snpy'):
    """
    Applys cuts to paramater file
    :param path: str; location of paramter file
    :param algo: str; algorithm to apply cuts
    """
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
def merge_snpy_salt_params(snpy_path: str, salt_path: str, save_loc: str):
    """
    Combines SNooPy and SALT params into a single file.
    :param snpy_path: str; location of snpy file
    :param salt_path: str; location of salt file
    :param save_loc: str; location of merged file
    """
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
def resid_v_z(path: str, title: str = '', save_loc: str = ''):
    """
    Plots the Hubble Residual v. Redshift
    :param path: str; location of file of parameter file
    :param title: str = ''; optional title to put at top of file
    :param save_loc: str = ''; location to save plot
    """
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
    if len(save_loc) != 0:
        print('Saved figure to... ', save_loc)
        plt.savefig(save_loc)
    plt.show()
    return
def resid_v_mass(path: str, title='', save_loc=None):
    """
    Plots the Hubble Residual v. Mass
    :param path: str; location of file of parameter file
    :param title: str = ''; optional title to put at top of file
    :param save_loc: str = ''; location to save plot
    """
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
def snpy_hist(hicat_params_file: str, norm_params_file: str, save_loc: str = '',
              sigma: float = 3.0, st_width: float = 0.02, c_width: float = 0.02):
    """
    Histogram of the SNooPy paramaters of 91bg-like vs. normal SNe Ia
    :param hicat_params_file: str; location of 91bg-like paramaters
    :param norm_params_file: str; location of normal paramaters
    :param save_loc: str = '';
    :param sigma: float = 3.0;
    :param st_width: float = 0.02;
    :param c_width: float = 0.02;
    """
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
def salt_hist(hicat_params_file: str, norm_params_file: str, save_loc: str = '',
              sigma: float = 3.0, st_width: float = 0.3, c_width: float = 0.08):
    """
    Histogram of the SALT paramaters of 91bg-like vs. normal SNe Ia
    :param hicat_params_file: str; location of 91bg-like paramaters
    :param norm_params_file: str; location of normal paramaters
    :param save_loc: str = '';
    :param sigma: float = 3.0;
    :param st_width: float = 0.02;
    :param c_width: float = 0.02;
    """
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
    """
    Updates the plots on the README.md file
    """
    # Update Mass Plots
    resid_v_mass(path='output/combiend__snpy_params_cut.txt',
                 title='Hubble Residual v. Host Stellar Mass of CSP-ATLAS-ZTF 91bg-like SNe Ia [SNooPy]',
                 save_loc='saved/readme_plots/csp-atlas-ztf_snpy_resid_v_mass.png')
    resid_v_mass(path='output/combiend__salt_params_cut.txt',
                 title='Hubble Residual v. Host Stellar Mass of CSP-ATLAS-ZTF 91bg-like SNe Ia [SALT3]',
                 save_loc='saved/readme_plots/csp-atlas-ztf_salt_resid_v_mass.png')
    resid_v_mass(path='output/merged_params_cut.txt',
                 title='Hubble Residual v. Host Stellar Mass of CSP-ATLAS-ZTF 91bg-like SNe Ia [SALT3-SNooPy]',
                 save_loc='saved/readme_plots/merged_resid_v_mass.png')
    resid_v_mass(path='txts/norm_10-11-24/normSNe_merged_10-11-24.txt',
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
def mass_step_calc(mu: numpy.array(float), mu_err: numpy.array(float), mass: numpy.array(float), z: numpy.array(float),
                   cut: float = 10.0) -> (dict, dict):
    """
    Calculates the mass step given arrays data
    :param mu: numpy.array(float);
    :param mu_err: numpy.array(float);
    :param mass: numpy.array(float);
    :param z: numpy.array(float);
    :param cut: float = 10.0;
    :return dict; two dictionary of mass step and error & lower/upper weighted averages
    """
    if cut == 'median':
        cut = round(np.median(mass), 4)

    resid = mu - gen.current_cosmo().distmod(z).value

    upper_resid = np.average(resid[mass > cut], weights=(1/(mu_err[mass > cut]**2)))
    lower_resid = np.average(resid[mass < cut], weights=(1/(mu_err[mass < cut]**2)))

    upper_resid_err = np.std(resid[mass > cut]) / np.sqrt(len(mu_err[mass > cut]))
    lower_resid_err = np.std(resid[mass < cut]) / np.sqrt(len(mu_err[mass < cut]))

    mass_step = np.abs(upper_resid - lower_resid)
    mass_step_err = np.sqrt((lower_resid_err**2) + (upper_resid_err**2))


    return ({'value': mass_step, 'err': mass_step_err},
            {'lower_resid': {'value': lower_resid, 'err': lower_resid_err},
             'upper_resid': {'value': upper_resid, 'err': upper_resid_err}})
def dataset_analysis():
    """
    Displays data to update Table 2 in Overleaf table
    """
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

# Main Function Call ------------------------------------------------------------------------------------------------ #
def main_help():
    """
    Call to get examples of calls for main()
    """
    print('===========================================================================================================')
    print('Ex. Individual:', "main(fit_type='indv', data_set='CSP', algo='snpy', path='data/CSP/SN2005ke_snpy.txt')")
    print('------')
    print('Ex. Batch:', "main(fit_type='batch', data_set='CSP', algo='snpy', dmag_max=1.00)")
    print('------')
    print('Ex. Combined:', "main(fit_type='combiend', algo='snpy', dmag_max=1.00)")
    print('===========================================================================================================')
    return
def main(fit_type: str, data_set: str = '', path: str = '', algo: str = '', save_loc: str = '',
         dmag_max: float = 0.00, dflux_max: float = 0.00) -> list[object]:
    """
    A function to easily fit data using both algorithms.
    :param fit_type: str; type of fitting protocol to use
    :param data_set = '': str; name of data set
    :param algo = '': str; algorithms to fit data with (snpy/salt)
    :param path = '': str; data location for individual fits (only call for fit_type='indv')
    :param save_loc = '': str; location to save parameter data of SNe fits
    :param dmag_max: float=0.00; magnitude error cut off for taking in data
    :param dflux_max: float=0.00; flux error cut off for taking in data
    :return: SNe: list[object]; List of sn91bg() classes from fitting call.
    """
    fit_type, data_set, algo = fit_type.lower(), data_set.lower(), algo.lower()  # Drops everything to lowercase
    SNe = []  # Initialize return array
    if fit_type == 'indv':
        SN = SN91bg(path=path, data_set = data_set, dmag_max = dmag_max, dflux_max = dflux_max)
        if len(SN.objname) != 0:
            SN.fit(algo)
            SNe.append(SN)
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

    main(fit_type='combiend', algo='snpy', dmag_max=1.00)
    main(fit_type='combiend', algo='salt', dmag_max=1.00)
    merge_snpy_salt_params('output/combiend__snpy_params_cut.txt',
                           'output/combiend__salt_params_cut.txt',
                           'output/merged_params_cut.txt')
    update_readme_plots()




    # # SNe = main(fit_type='combiend', algo='salt', dmag_max=1.00)
    # # save_params_to_file_cov('output/combiend__salt_params_cov.txt', SNe)
    # # sample_cutter('output/combiend__salt_params_cov.txt', 'salt')
    # opt_dict = optimize_alpha_beta('output/combiend__salt_params_cov_cut.txt')
    # print('91bg-like SNe Ia ========================')
    # print('\tAlpha: '+str(round(opt_dict['alpha'], 3)))
    # print('\tBeta:  '+str(round(opt_dict['beta'], 3)))
    #
    # # SNe = norm_fit(algo = 'salt', save_loc = 'txts/norm_10-18-24/norm_10-18-24_params.txt', dflux_max = 1.00)
    # # save_params_to_file_cov('txts/norm_10-18-24/norm_10-18-24_params_cov.txt', SNe)
    # # sample_cutter('txts/norm_10-18-24/norm_10-18-24_params_cov.txt', 'salt')
    # opt_dict = optimize_alpha_beta('txts/norm_10-18-24/norm_10-18-24_params_cov_cut.txt')
    # print('Normal SNe Ia ===========================')
    # print('\tAlpha: '+str(round(opt_dict['alpha'], 3)))
    # print('\tBeta:  '+str(round(opt_dict['beta'], 3)))
    #
    # print('Expected Normal SNe =====================')
    # print('\tAlpha: 0.153')
    # print('\tBeta:  2.980')

    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')
