# M.D. Woods - 10/17/24
import os
import sys
import numpy
import snpy
import glob
import corner
import shutil
import sncosmo
import numpy as np
import time as systime
import datetime
import matplotlib
import matplotlib.pyplot as plt
from astro_ghost.ghostHelperFunctions import getTransientHosts
from astropy.coordinates import SkyCoord, Galactic
from astropy.table import Table
from astroquery.sdss import SDSS
from astropy.stats import sigma_clip, sigma_clipped_stats

from scripts import general as gen
from scripts.salt3_param_fitter import optimize_alpha_beta
from scripts import get_vpec  # Loads heavy data so if not activly fitting I would comment this out

CONSTANTS = gen.get_constants()
CURRENTDATE = datetime.datetime.now()
COLOR_WHEEL = {'ZTF': '#D8973C', 'ATLAS': '#BD632F', 'CSP': '#72513B', 'ATLAS-ZTF': '#273E47',
               'ZTF_SNPY': '#D8973C', 'ATLAS_SNPY': '#BD632F', 'CSP_SNPY': '#72513B', 'ATLAS-ZTF_SNPY': '#273E47',
               'ZTF_SALT': '#D8973C', 'ATLAS_SALT': '#BD632F', 'CSP_SALT': '#72513B', 'ATLAS-ZTF_SALT': '#273E47',
               'Pan+': '#BD632F', 'PanPlus': '#BD632F', 'Histogram': '#3B5058',
               '10': '#A4243B', 'median': '#D8C99B'}

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
            self.origin, self.objname, self.originalname = hdr[0], hdr[1], hdr[2]
            self.coords, self.z, self.z_cmb = (float(hdr[3]), float(hdr[4])), float(hdr[5]), float(hdr[6])
            self.discovery_date = float(hdr[7][:-1])

            f.readline() # Skip +++++++ line

            # Read in parameter data
            line = f.readline()
            while '+++' not in line:
                line = line.split(', ')
                self.params.update({line[0]: {'value': float(line[1]), 'err': float(line[2])}})
                line = f.readline()

            f.readline() # Skip +++++++ line

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

        # Calculate average magnitude around peak
        window = 5
        if algo == "snpy":
            peak_time = self.params['Tmax']['value']
        elif algo == "salt":
            peak_time = self.params['t0']['value']
        self.params.update({'peak_mag': {'value': np.average(self.mag[(self.time < peak_time+window) & (self.time > peak_time-window)]),
                                         'err': np.average(self.dmag[(self.time < peak_time+window) & (self.time > peak_time-window)])}})

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

            mu = pho_mB + (alpha * self.params['x1']['value']) - (beta * self.params['c']['value']) - M0
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
        f.write(f'# Created by M.D. Woods -- {CURRENTDATE} -- NUM TARGETS: {len(SNe)}\n')
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
def sample_cutter(path: str, algo: str = 'snpy', save_loc: str = ''):
    """
    Applys cuts to paramater file
    :param path: str; location of paramter file
    :param algo: str; algorithm to apply cuts
    """
    hdr, data = gen.default_open(path)
    if len(data) == 0:
        print('[+++] No params to cut!')
        return
    original_num = str(len(data[:, 0]))
    if algo == 'snpy':
        print('[+++] Cutting sample for SNooPy data...')
        cuts = {'z': 0.015, 'EBVhost': (-0.2, 0.3), 'EBVhost_err': 0.1, 'st': (-999, 1.0), 'st_err': 0.1, 'Tmax_err': 1, 'mu_err': 0.2}
        # cuts = {'z': 0.015, 'EBVhost': (0.1, 999), 'EBVhost_err': 0.1, 'st': (-999, 1.0), 'st_err': 0.1, 'Tmax_err': 1, 'mu_err': 0.2}
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
        cuts = {'z': 0.015, 'c': (-0.6, 0.6), 'x1': (-4.5, 4.5), 'c_err': 0.1, 'x1_err': 1, 't0_err': 1, 'mu_err': 0.2}
        # cuts = {'z': 0.015, 'c': (0.1, 999), 'x1': (-4.5, 4.5), 'c_err': 0.1, 'x1_err': 1, 't0_err': 1, 'mu_err': 0.2}
        f_out = '      | '
        for c in cuts:
            f_out += c + ': ' + str(cuts[c]) + ' | '
        print(f_out)
        data = data[(data[:, hdr.index('z_cmb')].astype(float) > cuts['z']) &
                    (data[:, hdr.index('c')].astype(float) > cuts['c'][0]) &
                    (data[:, hdr.index('c')].astype(float) < cuts['c'][1]) &
                    (data[:, hdr.index('x1')].astype(float) > cuts['x1'][0]) &
                    (data[:, hdr.index('x1')].astype(float) < cuts['x1'][1]) &
                    (data[:, hdr.index('c_err')].astype(float) < cuts['c_err']) &
                    (data[:, hdr.index('x1_err')].astype(float) < cuts['x1_err']) &
                    (data[:, hdr.index('t0_err')].astype(float) < cuts['t0_err']) &
                    (data[:, hdr.index('mu_err')].astype(float) < cuts['mu_err'])]

    # # Remove mass outliers -- bad dont do this
    # mass = data[:, hdr.index('hostMass')].astype(float)
    # m_mn, m_md, m_std = sigma_clipped_stats(mass)
    # data = data[abs(mass - m_mn) < 3 * m_std]

    resid = (data[:, hdr.index('mu')].astype(float) -
             gen.current_cosmo().distmod(data[:, hdr.index('z_cmb')].astype(float)).value)
    r_mn, r_md, r_std = sigma_clipped_stats(resid)
    data = data[abs(resid - r_mn) < 3 * r_std]

    new_resid = (data[:, hdr.index('mu')].astype(float) -
             gen.current_cosmo().distmod(data[:, hdr.index('z_cmb')].astype(float)).value)
    mn, md, std = sigma_clipped_stats(new_resid)
    print('      Hubble Residual Scatter:', std)

    # Save to file
    if len(save_loc) == 0: save_loc = path[:-4]+'_cut.txt'
    with open(save_loc, 'w') as f:
        # Info
        f.write(f'# Created by M.D. Woods -- {CURRENTDATE} -- HUBBLE RESID STD: {std}\n')
        f.write(f'# {cuts}\n')

        # Header
        f_out = hdr[0]
        for h in hdr[1:]:
            f_out += ', ' + h
        f.write(f_out + '\n')

        # Data
        for line in data[:]:
            f_out = line[0]
            for item in line[1:]:
                f_out += ', ' + str(item)
            f.write(f_out + '\n')
    print(f'      Cut file saved to... {save_loc}')
    print(f'      [ {data.shape[0]} / {original_num} ]')

    return
def merge_params(snpy_path: str, salt_path: str, save_loc: str = '', algo_bias: str = 'best'):
    """
    Combines SNooPy and SALT params into a single file.
    :param snpy_path: str; location of snpy file
    :param salt_path: str; location of salt file
    :param save_loc: str; location of merged file
    :param algo_bias; str; algorithm to prioritize
    """
    # Load data
    hdr_snpy, snpy_data = gen.default_open(snpy_path)
    hdr_salt, salt_data = gen.default_open(salt_path)

    # Order the paramaters
    new_snpy_data = np.array([
        snpy_data[:, hdr_snpy.index('objname')],
        snpy_data[:, hdr_snpy.index('ra')],
        snpy_data[:, hdr_snpy.index('dec')],
        snpy_data[:, hdr_snpy.index('z')],
        snpy_data[:, hdr_snpy.index('z_cmb')],
        snpy_data[:, hdr_snpy.index('MJDs')],
        snpy_data[:, hdr_snpy.index('MJDe')],
        snpy_data[:, hdr_snpy.index('origin')],
        snpy_data[:, hdr_snpy.index('mu')],
        snpy_data[:, hdr_snpy.index('mu_err')],
        snpy_data[:, hdr_snpy.index('hostMass')],
        snpy_data[:, hdr_snpy.index('hostMass_err')],
        snpy_data[:, hdr_snpy.index('Tmax')],
        snpy_data[:, hdr_snpy.index('Tmax_err')],
        snpy_data[:, hdr_snpy.index('st')],
        snpy_data[:, hdr_snpy.index('st_err')],
        snpy_data[:, hdr_snpy.index('EBVhost')],
        snpy_data[:, hdr_snpy.index('EBVhost_err')],
        snpy_data[:, hdr_snpy.index('peak_mag')],
        snpy_data[:, hdr_snpy.index('peak_mag_err')]
    ])
    snpy_data = new_snpy_data.T
    new_salt_data = np.array([
        salt_data[:, hdr_salt.index('objname')],
        salt_data[:, hdr_salt.index('ra')],
        salt_data[:, hdr_salt.index('dec')],
        salt_data[:, hdr_salt.index('z')],
        salt_data[:, hdr_salt.index('z_cmb')],
        salt_data[:, hdr_salt.index('MJDs')],
        salt_data[:, hdr_salt.index('MJDe')],
        salt_data[:, hdr_salt.index('origin')],
        salt_data[:, hdr_salt.index('mu')],
        salt_data[:, hdr_salt.index('mu_err')],
        salt_data[:, hdr_salt.index('hostMass')],
        salt_data[:, hdr_salt.index('hostMass_err')],
        salt_data[:, hdr_salt.index('t0')],
        salt_data[:, hdr_salt.index('t0_err')],
        salt_data[:, hdr_salt.index('x1')],
        salt_data[:, hdr_salt.index('x1_err')],
        salt_data[:, hdr_salt.index('c')],
        salt_data[:, hdr_salt.index('c_err')],
        salt_data[:, hdr_salt.index('peak_mag')],
        salt_data[:, hdr_salt.index('peak_mag_err')]
    ])
    salt_data = new_salt_data.T

    # Tag arrays for sorting
    snpy_data = np.hstack((snpy_data, np.full((len(snpy_data[:, 0]), 1), 'SNPY')))  # Add SNooPy tag
    salt_data = np.hstack((salt_data, np.full((len(salt_data[:, 0]), 1), 'SALT')))  # Add SALT3 tag

    # Modify origins for plotting later
    snpy_data[:, 7] = np.char.add(snpy_data[:, 7],
                                  np.full(len(snpy_data[:, 7]), '_SNPY'))
    salt_data[:, 7] = np.char.add(salt_data[:, 7],
                                  np.full(len(salt_data[:, 7]), '_SALT'))

    # Make merged arrays
    combined_arr = np.vstack((snpy_data, salt_data))  # Stack arrays
    new_merged_arr = np.empty((1, len(combined_arr[0, :])))  # Make empty array for validated data

    # Sort through arrays
    print(f'Prioritizing for {algo_bias.upper()}...')
    for arr in combined_arr:
        if arr[0] in new_merged_arr[:, 0]:
            continue
        elif (algo_bias.upper() == 'BEST' and
                arr[0] in combined_arr[:, 0][combined_arr[:, -1] == 'SNPY'] and
                arr[0] in combined_arr[:, 0][combined_arr[:, -1] == 'SALT']):
            if (combined_arr[(combined_arr[:, -1] == 'SNPY') & (combined_arr[:, 0] == arr[0])][0, -4] <
                combined_arr[(combined_arr[:, -1] == 'SALT') & (combined_arr[:, 0] == arr[0])][0, -4]):
                new_merged_arr = np.vstack((new_merged_arr,
                                            combined_arr[(combined_arr[:, -1] == 'SNPY') & (combined_arr[:, 0] == arr[0])]))
            else:
                new_merged_arr = np.vstack((new_merged_arr,
                                            combined_arr[(combined_arr[:, -1] == 'SALT') & (combined_arr[:, 0] == arr[0])]))
        elif len(combined_arr[(combined_arr[:, -1] == algo_bias.upper()) & (combined_arr[:, 0] == arr[0])]) != 0:
            new_merged_arr = np.vstack((new_merged_arr,
                                        combined_arr[(combined_arr[:, -1] == algo_bias.upper()) & (combined_arr[:, 0] == arr[0])]))
        else:
            new_merged_arr = np.vstack((new_merged_arr, arr))
    new_merged_arr = new_merged_arr[1:, :] # Remove empty first array

    # Write to file
    if len(save_loc) != 0:
        with open(save_loc, 'w') as f:
            f.write(f'# Created by M.D. Woods -- {CURRENTDATE} -- NUM TARGETS: {len(new_merged_arr)}\n')
            f.write(f"# WARNING: 'stretch' and 'color' are NOT identical for SNPY & SALT.\n")
            print('objname, ra, dec, z, z_cmb, MJDs, MJDe, origin, mu, mu_err, hostMass, hostMass_err, Tmax, Tmax_err, '
                  'stretch, stretch_err, color, color_err, peak_mag, peak_mag_err, algo', file=f)
            for arr in new_merged_arr:
                print(f'{arr[0]}, {arr[1]}, {arr[2]}, {arr[3]}, {arr[4]}, {arr[5]}, {arr[6]}, {arr[7]}, {arr[8]}, '
                      f'{arr[9]}, {arr[10]}, {arr[11]},  {arr[12]}, {arr[13]}, {arr[14]}, {arr[15]},  {arr[16]}, '
                      f'{arr[17]}, {arr[18]}, {arr[19]}, {arr[-1]}', file=f)
        print('Saved merged param file to... ', save_loc)

    print(f'Merged SNooPy+SALT3 -- Total: {len(new_merged_arr[:, -1])} ['+
          f"SNooPy: {len(new_merged_arr[:, -1][new_merged_arr[:, -1] == 'SNPY'])}, "+
          f"SALT3: {len(new_merged_arr[:, -1][new_merged_arr[:, -1] == 'SALT'])}]")
    return new_merged_arr
    # return
def merged_options(path_snpy: str, path_salt: str, save_loc: str = '', plot: bool = False):
    all_std, all_choices = [], ['snpy', 'salt', 'best']
    for choice in all_choices:
        data = merge_params(path_snpy,
                            path_salt,
                            algo_bias=choice)
        mn, md, std = sigma_clipped_stats(gen.get_resid(data[:, 8], data[:, 4]))
        all_std.append(std)

        if plot:
            resid_v_mass(save_loc)

    print(f'[+++] Lowest Residual STD Option: {all_choices[all_std.index(min(all_std))]}')
    merge_params(path_snpy,
                 path_salt,
                 save_loc,
                 algo_bias=all_choices[all_std.index(min(all_std))])

    return

# Plotting Functions ------------------------------------------------------------------------------------------------ #
def resid_v_z(path: str, title: str = '', save_loc: str = ''):
    """
    Plots the Hubble Residual v. Redshift
    :param path: str; location of file of parameter file
    :param title: str = ''; optional title to put at top of file
    :param save_loc: str = ''; location to save plot
    """
    fig, axs = plt.subplots(1, 2, figsize=(15, 6), gridspec_kw={'width_ratios': [10, 1]}, constrained_layout=True)

    # Pull data from saved text & header
    hdr, data = gen.default_open(path)

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
                        color=COLOR_WHEEL[origin], elinewidth=0.8, **format_dict)

    # Make histogram
    axs[1].hist(resid_mu, bins=int((np.max(resid_mu) - np.min(resid_mu)) / 0.1),  # Bin Width = 0.1
                orientation="horizontal", color=COLOR_WHEEL['Histogram'])

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
def resid_v_mass(path: str, title: str = '', cuts: list = [10, 'median'], save_loc: str = '', label: bool = False):
    """
    Plots the Hubble Residual v. Mass
    :param path: str; location of file of parameter file
    :param title: str = ''; optional title to put at top of file
    :param save_loc: str = ''; location to save plot
    """
    fig, axs = plt.subplots(1, 2, figsize=(15, 6), gridspec_kw={'width_ratios': [10, 1]}, constrained_layout=True)

    # Pull data from saved text & header
    hdr, data = gen.default_open(path)

    # Set Arrays
    names = data[:, hdr.index('objname')]
    z = data[:, hdr.index('z_cmb')].astype(float)
    origins = data[:, hdr.index('origin')]
    mass = data[:, hdr.index('hostMass')].astype(float)
    mass_err =  data[:, hdr.index('hostMass_err')].astype(float)
    mu = data[:, hdr.index('mu')].astype(float)
    mu_err = data[:, hdr.index('mu_err')].astype(float)

    # Calculate Hubble Residual
    resid_mu = mu - gen.current_cosmo().distmod(z).value
    resid_mu_err = np.copy(mu_err)

    # Subtracting off Average Hubble Residual
    resid_mu -= np.average(resid_mu[~np.isnan(resid_mu)])

    # Corrections
    mu_err = np.sqrt(mu_err ** 2.0 + 0.1 ** 2.0)  # intrinsic dispersion added in quadrature
    # r_mn, r_md, r_std = sigma_clipped_stats(resid_mu)
    # m_mn, m_md, m_std = sigma_clipped_stats(mass)
    # iGood = (abs(resid_mu - r_mn) < 3 * r_std) & (abs(mass - m_mn) < 3 * m_std)
    # mu, mu_err, z, resid_mu, resid_mu_err, mass, mass_err, origins, names = (
    #     mu[iGood], mu_err[iGood], z[iGood], resid_mu[iGood], resid_mu_err[iGood], mass[iGood], mass_err[iGood], origins[iGood], names[iGood])

    # Make main plot
    for origin in np.unique(origins):
        format_dict = {'marker': 'o', 'fmt': 'o', 'label': origin, 'alpha': 1, 'ms': 6}
        if 'algo' in hdr:
            if 'SNPY' in origin:
                format_dict['label'], format_dict['marker'] = origin[:-5], 'o'
            elif 'SALT' in origin:
                format_dict['label'], format_dict['marker'] = origin[:-5], '^'

        indexes = np.where(origins == origin)[0]
        axs[0].errorbar(mass[indexes], resid_mu[indexes], xerr=mass_err[indexes], yerr=resid_mu_err[indexes],
                        color=COLOR_WHEEL[origin], elinewidth=0.8, **format_dict)

    # Labels
    if label:
        for i in range(len(resid_mu)):
            axs[0].text(mass[i], resid_mu[i], names[i], ha='left', va='top', size='xx-small')

    # Make histogram
    axs[1].hist(resid_mu, bins=int((np.max(resid_mu) - np.min(resid_mu)) / 0.05),  # Bin Width = 0.05
                orientation="horizontal", color=COLOR_WHEEL['Histogram'])

    # Extra Info
    extra_info = '$\sigma$: '+str(round(np.std(resid_mu), 4)) + ', $n$: ' + str(len(resid_mu))
    if 'merged' in path:
        extra_info += r' | SALT3: $\triangle$, SNooPy: $\bigcirc$'
    axs[0].text(np.min(mass), np.max(resid_mu),
                extra_info,
                horizontalalignment='left', verticalalignment='bottom')

    # Display Both Mass Steps
    for cut in cuts:
        num_cut, lin_color = cut, COLOR_WHEEL[str(cut)]
        if cut == 'median':
            num_cut = round(np.median(mass), 4)

        # Get Mass Step
        mass_step_dict, resid_dict = mass_step_calc(mu, mu_err, resid_mu, mass, z, cut=num_cut)
        if mass_step_dict['value'] == 0.00 and mass_step_dict['err'] == 0.00:
            continue

        # Plot Mass Step
        lin_details = {'linestyle': '--', 'linewidth': 1.0, 'color': lin_color}
        fill_details = {'color': lin_color, 'alpha': 0.15}
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

    # Labels
    fig.suptitle(title)
    axs[0].set(xlabel="Host Stellar Mass ($\log M_{*}[M_{\odot}]$)",
               ylabel='Hubble Residuals (mag)')  # Sub-plot Labels
    axs[1].get_yaxis().set_visible(False)  # Turn off y-axis labels
    axs[0].legend(loc='best')

    # Adjust axies
    clearance = 0.2
    axs[0].set_ylim(np.min(resid_mu) - clearance, np.max(resid_mu) + clearance)
    axs[1].set_ylim(np.min(resid_mu) - clearance, np.max(resid_mu) + clearance)

    # Saving Figure
    if len(save_loc) != 0:
        print('Saved figure to... ', save_loc)
        plt.savefig(save_loc)
    plt.show()
    return
def resid_v_mass_colormap(path: str, param: str, title: str = '', cuts: list = [10, 'median'], save_loc: str = '',
                          label: bool = False):
    fig, axs = plt.subplots(1, 2, figsize=(15, 6), gridspec_kw={'width_ratios': [10, 1]},
                            constrained_layout=True)

    # Pull data from saved text & header
    hdr, data = gen.default_open(path)

    # Set Arrays
    names = data[:, hdr.index('objname')]
    z = data[:, hdr.index('z_cmb')].astype(float)
    origins = data[:, hdr.index('origin')]
    mass = data[:, hdr.index('hostMass')].astype(float)
    mass_err =  data[:, hdr.index('hostMass_err')].astype(float)
    mu = data[:, hdr.index('mu')].astype(float)
    mu_err = data[:, hdr.index('mu_err')].astype(float)

    # Select param
    try:
        param_arr = data[:, hdr.index(param)].astype(float)
    except IndexError:
        print(f"[!!!] '{param}' is not a valid paramater for '{path}'")

    # Calculate Hubble Residual
    resid_mu = mu - gen.current_cosmo().distmod(z).value
    resid_mu_err = np.copy(mu_err)

    # Subtracting off Average Hubble Residual
    resid_mu -= np.average(resid_mu[~np.isnan(resid_mu)])

    # Corrections
    mu_err = np.sqrt(mu_err ** 2.0 + 0.1 ** 2.0)  # intrinsic dispersion added in quadrature

    # Make main plot
    axs[0].errorbar(mass, resid_mu, xerr=mass_err, yerr=resid_mu_err, fmt='.')
    mass_plt = axs[0].scatter(mass, resid_mu, c=param_arr, cmap='viridis_r', label=f"CMAP: '{param}'")

    # Labels
    if label:
        for i in range(len(resid_mu)):
            axs[0].text(mass[i], resid_mu[i], names[i], ha='left', va='top', size='xx-small')

    # Make histogram
    axs[1].hist(resid_mu, bins=int((np.max(resid_mu) - np.min(resid_mu)) / 0.05),  # Bin Width = 0.05
                orientation="horizontal", color=COLOR_WHEEL['Histogram'])

    # Extra Info
    extra_info = '$\sigma$: '+str(round(np.std(resid_mu), 4)) + ', $n$: ' + str(len(resid_mu))
    if 'merged' in path:
        extra_info += r' | SALT3: $\triangle$, SNooPy: $\bigcirc$'
    axs[0].text(np.min(mass), np.max(resid_mu),
                extra_info,
                horizontalalignment='left', verticalalignment='bottom')

    # Display Both Mass Steps
    for cut in cuts:
        num_cut, lin_color = cut, COLOR_WHEEL[str(cut)]
        if cut == 'median':
            num_cut = round(np.median(mass), 4)

        # Get Mass Step
        mass_step_dict, resid_dict = mass_step_calc(mu, mu_err, resid_mu, mass, z, cut=num_cut)
        if mass_step_dict['value'] == 0.00 and mass_step_dict['err'] == 0.00:
            continue

        # Plot Mass Step
        lin_details = {'linestyle': '--', 'linewidth': 1.0, 'color': lin_color}
        fill_details = {'color': lin_color, 'alpha': 0.15}
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

    # Labels
    fig.suptitle(title)
    axs[0].set(xlabel="Host Stellar Mass ($\log M_{*}[M_{\odot}]$)",
               ylabel='Hubble Residuals (mag)')  # Sub-plot Labels
    axs[1].get_yaxis().set_visible(False)  # Turn off y-axis labels
    axs[0].legend(loc='best')
    # plt.colorbar(mass_plt, ax=axs[1])

    # Adjust axies
    clearance = 0.2
    axs[0].set_ylim(np.min(resid_mu) - clearance, np.max(resid_mu) + clearance)
    axs[1].set_ylim(np.min(resid_mu) - clearance, np.max(resid_mu) + clearance)

    # Saving Figure
    if len(save_loc) != 0:
        print('Saved figure to... ', save_loc)
        plt.savefig(save_loc)
    plt.show()
    return
def param_hist(hicat_params_file: str, norm_params_file: str, algo: str, type: str, save_loc: str = '',
    sigma: float = 3.0, st_width: float = 0.04, c_width: float = 0.04):
    """
    Histogram of the SNooPy paramaters of 91bg-like vs. normal SNe Ia
    :param hicat_params_file: str; location of 91bg-like paramaters
    :param norm_params_file: str; location of normal paramaters
    :param algo: str; algo of histogram to generate
    :param save_loc: str = '';
    :param type: str = 'median' or 'average';
    :param sigma: float = 3.0;
    :param st_width: float = 0.02;
    :param c_width: float = 0.02;
    """
    # Open data
    hicat_hdr, hicat_data = gen.default_open(hicat_params_file)
    norm_hdr, norm_data = gen.default_open(norm_params_file)

    # Select Algo
    if algo == 'snpy':
        param_dict = {'stretch': 'st', 'color': 'EBVhost'}
    elif algo == 'salt':
        param_dict = {'stretch': 'x1', 'color': 'c'}
    else:
        raise ValueError(f"[!!!] '{algo}' is invalid; choose between 'snpy' or 'salt'")

    fig, ax = plt.subplots(1, 2, figsize=(16, 4), constrained_layout=True)

    # Stretch plot
    st_norm = sigma_clip(norm_data[:, norm_hdr.index(param_dict['stretch'])].astype(float), sigma=sigma)
    st_hicat = sigma_clip(hicat_data[:, hicat_hdr.index(param_dict['stretch'])].astype(float), sigma=sigma)
    st_norm_err = sigma_clip(norm_data[:, norm_hdr.index(param_dict['stretch']+'_err')].astype(float), sigma=sigma)
    st_hicat_err = sigma_clip(hicat_data[:, hicat_hdr.index(param_dict['stretch']+'_err')].astype(float), sigma=sigma)
    ax[0].hist(st_norm, int((np.max(st_norm) - np.min(st_norm)) / st_width), label='Normal', color='#5AD2F4')
    ax[0].hist(st_hicat, int((np.max(st_hicat) - np.min(st_hicat)) / st_width), label='91bg', color='#62BEC1', alpha=0.75)
    ax[0].set_title(f"Stretch [{param_dict['stretch']}]")

    # Color plot
    c_norm = sigma_clip(norm_data[:, norm_hdr.index(param_dict['color'])].astype(float), sigma=sigma)
    c_hicat = sigma_clip(hicat_data[:, hicat_hdr.index(param_dict['color'])].astype(float), sigma=sigma)
    c_norm_err = sigma_clip(norm_data[:, norm_hdr.index(param_dict['color']+'_err')].astype(float), sigma=sigma)
    c_hicat_err = sigma_clip(hicat_data[:, hicat_hdr.index(param_dict['color']+'_err')].astype(float), sigma=sigma)
    ax[1].hist(c_norm, int((np.max(c_norm) - np.min(c_norm)) / c_width), label='Normal', color='#5AD2F4')
    ax[1].hist(c_hicat, int((np.max(c_hicat) - np.min(c_hicat)) / c_width), label='91bg', color='#62BEC1', alpha=0.75)
    ax[1].set_title(f"Color [{param_dict['color']}]")

    # Plot labels
    if type == 'median':
        ax[0].axvline(x=np.median(st_norm),
                      label=r'$\tilde{x}_{Normal}$ = ' + str(round(np.median(st_norm), 2))+' $\pm$ '+str(round(np.median(st_norm_err), 3)),
                      linewidth=2.5, color='#4bb0cc', linestyle='--')
        ax[0].axvline(x=np.median(st_hicat),
                      label=r'$\tilde{x}_{91bg}$ = ' + str(round(np.median(st_hicat), 2))+' $\pm$ '+str(round(np.median(st_hicat_err), 3)),
                      linewidth=2.5, color='#52a1a3', linestyle=':')
        ax[1].axvline(x=np.median(c_norm),
                      label=r'$\tilde{x}_{Normal}$ = ' + str(round(np.median(c_norm), 2))+' $\pm$ '+str(round(np.median(c_norm_err), 3)),
                      linewidth=2.5, color='#4bb0cc', linestyle='--')
        ax[1].axvline(x=np.median(c_hicat),
                      label=r'$\tilde{x}_{91bg}$ = '+str(round(np.median(c_hicat), 2))+' $\pm$ '+str(round(np.median(c_hicat_err), 3)),
                      linewidth=2.5, color='#52a1a3', linestyle=':')
    elif type == 'average':
        ax[0].axvline(x=np.average(st_norm),
                      label=r'$\tilde{x}_{Normal}$ = ' + str(round(np.average(st_norm), 2))+' $\pm$ '+str(round(np.average(st_norm_err), 3)),
                      linewidth=2.5, color='#4bb0cc', linestyle='--')
        ax[0].axvline(x=np.average(st_hicat),
                      label=r'$\tilde{x}_{91bg}$ = ' + str(round(np.average(st_hicat), 2))+' $\pm$ '+str(round(np.average(st_hicat_err), 3)),
                      linewidth=2.5, color='#52a1a3', linestyle=':')
        ax[1].axvline(x=np.average(c_norm),
                      label=r'$\tilde{x}_{Normal}$ = ' + str(round(np.average(c_norm), 2))+' $\pm$ '+str(round(np.average(c_norm_err), 3)),
                      linewidth=2.5, color='#4bb0cc', linestyle='--')
        ax[1].axvline(x=np.average(c_hicat),
                      label=r'$\tilde{x}_{91bg}$ = '+str(round(np.average(c_hicat), 2))+' $\pm$ '+str(round(np.average(c_hicat_err), 3)),
                      linewidth=2.5, color='#52a1a3', linestyle=':')
    else:
        raise ValueError(f"[!!!] '{type}' is invalid; choose between 'median' or 'average'")
    ax[0].legend()
    ax[1].legend()


    if len(save_loc) > 0:
        plt.savefig(save_loc, dpi=300)
    plt.show()
    return
def param_corner_plot(snpy_path: str, salt_path: str, save_loc: str = ''):
    data = merge_params(snpy_path, salt_path)

    # Plot all common params
    all_params = np.array([data[:, 8].astype(float), data[:, 10].astype(float), data[:, 12].astype(float)])
    figure = corner.corner(all_params.T,
                           labels=[r"$\mu$", r"$\log M/M_{\ast}$", r"$T_{Max}$"],
                           quantiles=[0.16, 0.5, 0.84],
                           show_titles=True,
                           title_kwargs={"fontsize": 12})
    if len(save_loc) != 0:
        plt.savefig(save_loc+'allparam_cornerPlot.png')
    plt.show()

    # Plot algo-specific params
    for algo in ['SNPY', 'SALT']:
        all_params = np.array([data[data[:, -1] == algo][:, 8].astype(float),
                               data[data[:, -1] == algo][:, 10].astype(float),
                               data[data[:, -1] == algo][:, 12].astype(float),
                               data[data[:, -1] == algo][:, 14].astype(float),
                               data[data[:, -1] == algo][:, 16].astype(float)])
        figure = corner.corner(all_params.T,
                               labels=[r"$\mu$",
                                       r"$\log M/M_{\ast}$",
                                       r"$T_{Max}$",
                                       r"Stretch",
                                       r"Color"],
                               quantiles=[0.16, 0.5, 0.84],
                               show_titles=True,
                               title_kwargs={"fontsize": 12})
        if len(save_loc) != 0:
            plt.savefig(save_loc+algo.lower()+'_cornerPlot.png')
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
def mass_step_calc(mu: numpy.array(float), mu_err: numpy.array(float), resid: numpy.array(float),
                   mass: numpy.array(float), z: numpy.array(float),
                   cut: float = 10.0) -> (dict, dict):
    """
    Calculates the mass step given arrays data
    :param mu: numpy.array(float);
    :param mu_err: numpy.array(float);
    :param resid: numpy.array(float);
    :param mass: numpy.array(float);
    :param z: numpy.array(float);
    :param cut: float = 10.0;
    :return dict; two dictionary of mass step and error & lower/upper weighted averages
    """
    if cut == 'median':
        cut = round(np.median(mass), 4)

    try:
        upper_resid = np.average(resid[mass > cut], weights=(1/(mu_err[mass > cut]**2)))
        lower_resid = np.average(resid[mass < cut], weights=(1/(mu_err[mass < cut]**2)))

        upper_resid_err = np.std(resid[mass > cut]) / np.sqrt(len(mu_err[mass > cut]))
        lower_resid_err = np.std(resid[mass < cut]) / np.sqrt(len(mu_err[mass < cut]))

        mass_step = np.abs(upper_resid - lower_resid)
        mass_step_err = np.sqrt((lower_resid_err**2) + (upper_resid_err**2))
    except ZeroDivisionError:
        return ({'value': 0.00, 'err': 0.00},
                {'lower_resid': {'value': 0.00, 'err': 0.00},
                 'upper_resid': {'value': 0.00, 'err': 0.00}})


    return ({'value': mass_step, 'err': mass_step_err},
            {'lower_resid': {'value': lower_resid, 'err': lower_resid_err},
             'upper_resid': {'value': upper_resid, 'err': upper_resid_err}})
def SNooPy_SALT3_overlap(snpy_path: str, salt_path: str):
    """
    Compares SNooPy and SALT3 overlap
    :param snpy_path:
    :param salt_path:
    """
    snpy_data = np.genfromtxt(snpy_path, delimiter=', ', skip_header=1, dtype=str)
    salt_data = np.genfromtxt(salt_path, delimiter=', ', skip_header=1, dtype=str)
    with open(snpy_path, 'r') as f:
        hdr_snpy = f.readline().split(', ')
        hdr_snpy[-1] = hdr_snpy[-1][:-1]
    with open(salt_path, 'r') as f:
        hdr_salt = f.readline().split(', ')
        hdr_salt[-1] = hdr_salt[-1][:-1]

    overlap = []
    for n in snpy_data[:, 0]:
        if n in salt_data[:, 0]:
            overlap.append(n)

    for n_ind in range(len(overlap)):
        print('--------------\n', overlap[n_ind], '\n--------------')
        snpy_resid_mu  = (float(snpy_data[snpy_data[:, 0] == overlap[n_ind]][0, hdr_snpy.index('mu')]) -
                           gen.current_cosmo().distmod(float(snpy_data[snpy_data[:, 0] == overlap[n_ind]][0, hdr_snpy.index('z_cmb')])).value)
        snpy_resid_mu_err = float(snpy_data[snpy_data[:, 0] == overlap[n_ind]][0, hdr_snpy.index('mu_err')])
        print(f"{round(snpy_resid_mu,3)} +/- {round(snpy_resid_mu_err,3)}")

        salt_resid_mu = (float(salt_data[salt_data[:, 0] == overlap[n_ind]][0, hdr_salt.index('mu')]) -
                         gen.current_cosmo().distmod(float(salt_data[salt_data[:, 0] == overlap[n_ind]][0, hdr_salt.index('z_cmb')])).value)
        salt_resid_mu_err = float(salt_data[salt_data[:, 0] == overlap[n_ind]][0, hdr_salt.index('mu_err')])
        print(f"{round(salt_resid_mu, 3)} +/- {round(salt_resid_mu_err, 3)}")

        plt.errorbar(snpy_resid_mu, salt_resid_mu, xerr=snpy_resid_mu_err, yerr=salt_resid_mu_err, fmt='o')
    plt.xlabel('SNooPy'); plt.ylabel('SAlT3')
    plt.show()

    return
def alpha_beta_fitting():
    # SNe = main(fit_type='combiend', algo='salt', dmag_max=1.00)
    # save_params_to_file_cov('output/combiend__salt_params_cov.txt', SNe)
    # sample_cutter('output/combiend__salt_params_cov.txt', 'salt')
    opt_dict = optimize_alpha_beta('output/old/combiend__salt_params_cov_cut.txt')
    print('91bg-like SNe Ia ========================')
    print('\tAlpha: '+str(round(opt_dict['alpha'], 3)))
    print('\tBeta:  '+str(round(opt_dict['beta'], 3)))

    # SNe = norm_fit(algo = 'salt', save_loc = 'txts/norm_10-18-24/norm_10-18-24_params.txt', dflux_max = 1.00)
    # save_params_to_file_cov('txts/norm_10-18-24/norm_10-18-24_params_cov.txt', SNe)
    # sample_cutter('txts/norm_10-18-24/norm_10-18-24_params_cov.txt', 'salt')
    opt_dict = optimize_alpha_beta('txts/norm_10-18-24/norm_10-18-24_params_cov_cut.txt')
    print('Normal SNe Ia ===========================')
    print('\tAlpha: '+str(round(opt_dict['alpha'], 3)))
    print('\tBeta:  '+str(round(opt_dict['beta'], 3)))

    print('Expected Normal SNe =====================')
    print('\tAlpha: 0.153')
    print('\tBeta:  2.980')
    return
def data_stats(fmt: str = 'print', path: str = 'output/merged_params_cut.txt'):
    hdr, data = gen.default_open(path)

    all_origins = data[:, hdr.index('origin')]
    for i in range(len(all_origins)):
        all_origins[i] = all_origins[i][:-5]

    if fmt == 'print':
        # Print out format
        print(f"Origin, Number of SNe, CMB Redshift Range, Decl. Range, Avg. Mag. (+/- 5MJD), Avg. Mag. Err (+/- 5MJD)")
        for origin in np.unique(all_origins):
            indexs = np.where(all_origins == origin)[0]
            local_data = data[indexs]
            z = local_data[:, hdr.index('z_cmb')].astype(float)
            dec = local_data[:, hdr.index('dec')].astype(float)
            peak_mag = local_data[:, hdr.index('peak_mag')].astype(float)
            peak_mag_err = local_data[:, hdr.index('peak_mag_err')].astype(float)
            print(f"{origin}, "
                  f"{len(indexs)}, "
                  f"{round(np.min(z), 4)} to {round(np.max(z), 4)}, "
                  f"{round(np.min(dec), 4)} to {round(np.max(dec), 4)}, "
                  f"{round(np.average(peak_mag[~np.isnan(peak_mag)]), 4)}, "
                  f"{round(np.average(peak_mag_err[~np.isnan(peak_mag_err)]), 4)}")
        print(f"All Surveys,  "
              f"{len(data[:, 0])}, "
              f"{round(np.min(data[:, hdr.index('z_cmb')].astype(float)), 4)} to {round(np.max(data[:, hdr.index('z_cmb')].astype(float)), 4)}, "
              f"{round(np.min(data[:, hdr.index('dec')].astype(float)), 4)} to {round(np.max(data[:, hdr.index('dec')].astype(float)), 4)}, "
              f"{round(np.average(data[:, hdr.index('peak_mag')].astype(float)[~np.isnan(data[:, hdr.index('peak_mag')].astype(float))]), 4)}, "
              f"{round(np.average(data[:, hdr.index('peak_mag_err')].astype(float)[~np.isnan(data[:, hdr.index('peak_mag_err')].astype(float))]), 4)}")
    elif fmt == 'latex':
        # LaTeX format
        for origin in np.unique(all_origins):
            indexs = np.where(all_origins == origin)[0]
            local_data = data[indexs]
            z = local_data[:, hdr.index('z_cmb')].astype(float)
            dec = local_data[:, hdr.index('dec')].astype(float)
            peak_mag = local_data[:, hdr.index('peak_mag')].astype(float)
            peak_mag_err = local_data[:, hdr.index('peak_mag_err')].astype(float)
            print(f"{origin} & "
                  f"${len(indexs)}$ & "
                  f"${round(np.min(z), 4)}$ & ${round(np.max(z), 4)}$ & "
                  f"${round(np.min(dec), 4)}$ & ${round(np.max(dec), 4)}$ & "
                  f"${round(np.average(peak_mag[~np.isnan(peak_mag)]), 4)}$ & "
                  f"${round(np.average(peak_mag_err[~np.isnan(peak_mag_err)]), 4)}$\\\\")
        print(f"All Surveys &  "
              f"${len(data[:, 0])}$ & "
              f"${round(np.min(data[:, hdr.index('z_cmb')].astype(float)), 4)}$ & ${round(np.max(data[:, hdr.index('z_cmb')].astype(float)), 4)}$ & "
              f"${round(np.min(data[:, hdr.index('dec')].astype(float)), 4)}$ & ${round(np.max(data[:, hdr.index('dec')].astype(float)), 4)}$ & "
              f"${round(np.average(data[:, hdr.index('peak_mag')].astype(float)[~np.isnan(data[:, hdr.index('peak_mag')].astype(float))]), 4)}$ & "
              f"${round(np.average(data[:, hdr.index('peak_mag_err')].astype(float)[~np.isnan(data[:, hdr.index('peak_mag_err')].astype(float))]), 4)}$\\\\")
    else:
        raise ValueError(f"[!!!] Unrecognized format: {fmt}, select 'print' or 'latex'")

    return

# Reformatting Functions -------------------------------------------------------------------------------------------- #
def format_panthplus(path: str = 'txts/PanPlus_Latest.FITRES', save_loc: str = 'output/panthplus_params.txt'):
    # Load data from PanPlus
    hdr = ('CID CIDint IDSURVEY TYPE FIELD CUTFLAG_SNANA zHEL zHELERR zCMB zCMBERR zHD zHDERR VPEC VPECERR MWEBV '
           'HOST_LOGMASS HOST_LOGMASS_ERR HOST_sSFR HOST_sSFR_ERR PKMJDINI SNRMAX1 SNRMAX2 SNRMAX3 PKMJD PKMJDERR '
           'x1 x1ERR c cERR mB mBERR x0 x0ERR COV_x1_c COV_x1_x0 COV_c_x0 NDOF FITCHI2 FITPROB RA DEC HOST_RA HOST_DEC '
           'HOST_ANGSEP TGAPMAX TrestMIN TrestMAX ELU HOSTGAL_SFR HOSTGAL_SFR_ERR HOSTGAL_sSFR HOSTGAL_sSFR_ERR CUTMASK '
           'MU MUMODEL MUERR MUERR_RENORM MUERR_RAW MUERR_VPEC MURES MUPULL M0DIF M0DIFERR CHI2 biasCor_mu biasCorErr_mu '
           'biasCor_mB biasCor_x1 biasCor_c biasScale_muCOV IDSAMPLE IZBIN').split(' ')
    data = np.genfromtxt(path, dtype=str)
    data = data[:, 1:] # Removes 'SN:' from rows

    # Corrections
    data = data[np.abs(data[:, hdr.index('HOST_LOGMASS_ERR')].astype(float)) < 1]
    data = data[data[:, hdr.index('HOST_LOGMASS')].astype(float) > 0]

    # Write to new file
    with open(save_loc, 'w') as f:
        f.write(f'# Created by M.D. Woods -- {CURRENTDATE} -- NUM TARGETS: {len(data[:, 0])}\n')
        f.write(f'objname, ra, dec, z, z_cmb, MJDs, MJDe, origin, mu, mu_err, x1, x1_err, t0, t0_err, c, c_err, chisquare, chisquare_err, hostMass, hostMass_err\n')
        for i in range(len(data[:, 0])):
            f.write(f"{data[i, hdr.index('CID')]}, "
                    f"{data[i, hdr.index('RA')]}, "
                    f"{data[i, hdr.index('DEC')]}, "
                    f"{data[i, hdr.index('zHEL')]}, "
                    f"{data[i, hdr.index('zCMB')]}, "
                    f"{np.nan}, "
                    f"{np.nan}, "
                    f"PanPlus, "
                    f"{data[i, hdr.index('MU')]}, "
                    f"{abs(data[i, hdr.index('MUERR')].astype(float))}, "
                    f"{data[i, hdr.index('x1')]}, "
                    f"{abs(data[i, hdr.index('x1ERR')].astype(float))}, "
                    f"{data[i, hdr.index('PKMJD')]}, "
                    f"{abs(data[i, hdr.index('PKMJDERR')].astype(float))}, "
                    f"{data[i, hdr.index('c')]}, "
                    f"{abs(data[i, hdr.index('cERR')].astype(float))}, "
                    f"{data[i, hdr.index('CHI2')]}, "
                    f"{np.nan}, "
                    f"{data[i, hdr.index('HOST_LOGMASS')]}, "
                    f"{abs(data[i, hdr.index('HOST_LOGMASS_ERR')].astype(float))}\n")
    return
def format_dr3(path: str = 'txts/DR3_fits.dat', save_loc: str = 'output/dr3_params.txt'):
    # Load Data
    data = np.genfromtxt(path, dtype=str)
    hdr, data = list(data[0, :]), data[1:, :]

    # Write to file
    with open(save_loc, 'w') as f:
        f.write(f'# Created by M.D. Woods -- {CURRENTDATE} -- NUM TARGETS: {len(data[:, 0])}\n')
        f.write(f'objname, z, z_cmb, origin, st, st_err, Tmax, Tmax_err, EBVhost, EBVhost_err\n')
        for i in range(len(data[:, 0])):
            f.write(f"{data[i, hdr.index('SN')]}, "
                    f"{data[i, hdr.index('zhelio')]}, "
                    f"{data[i, hdr.index('zcmb')]}, "
                    f"DR3, "
                    f"{data[i, hdr.index('st')]}, "
                    f"{abs(data[i, hdr.index('e_st')].astype(float))}, "
                    f"{data[i, hdr.index('Tmax')]}, "
                    f"{abs(data[i, hdr.index('eTmax')].astype(float))}, "
                    f"{data[i, hdr.index('EBVhost')]}, "
                    f"{abs(data[i, hdr.index('e_EBVhost')].astype(float))}\n")
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
def refresh_all(new_data: bool = False):
    """
    Updates the plots on the README.md file
    """
    if new_data:
        main(fit_type='combiend', algo='snpy', dmag_max=1.00)
        main(fit_type='combiend', algo='salt', dmag_max=1.00)
        norm_fit('snpy', 'output/norm_snpy_params.txt', dmag_max=1.00)
        norm_fit('salt', 'output/norm_salt_params.txt', dmag_max=1.00)

    # Merge raw files
    merged_options('output/combiend__snpy_params.txt',
                   'output/combiend__salt_params.txt',
                   'output/merged_params.txt')
    merged_options('output/norm_snpy_params.txt',
                   'output/norm_salt_params.txt',
                   'output/norm_merged_params.txt')

    # Make sure most recent cuts are applied
    sample_cutter('output/combiend__snpy_params.txt', 'snpy')
    sample_cutter('output/combiend__salt_params.txt', 'salt')
    merged_options('output/combiend__snpy_params_cut.txt',
                   'output/combiend__salt_params_cut.txt',
                   'output/merged_params_cut.txt')
    sample_cutter('output/norm_snpy_params.txt', 'snpy')
    sample_cutter('output/norm_salt_params.txt', 'salt')
    merged_options('output/norm_snpy_params.txt',
                   'output/norm_salt_params.txt',
                   'output/norm_merged_params_cut.txt')

    # Update Mass Plots
    resid_v_mass(path='output/combiend__snpy_params_cut.txt',
                 save_loc = 'saved/readme_plots/csp-atlas-ztf_snpy_resid_v_mass.png')
                 # title='Hubble Residual v. Host Stellar Mass of CSP-ATLAS-ZTF 91bg-like SNe Ia [SNooPy]',
    resid_v_mass(path='output/combiend__salt_params_cut.txt',
                 save_loc='saved/readme_plots/csp-atlas-ztf_salt_resid_v_mass.png')
                 # title='Hubble Residual v. Host Stellar Mass of CSP-ATLAS-ZTF 91bg-like SNe Ia [SALT3]',)
    resid_v_mass(path='output/merged_params_cut.txt',
                 save_loc='saved/readme_plots/merged_resid_v_mass.png')
                 # title='Hubble Residual v. Host Stellar Mass of CSP-ATLAS-ZTF 91bg-like SNe Ia [SALT3-SNooPy]')
    resid_v_mass(path='output/norm_merged_params_cut.txt',
                 save_loc='saved/readme_plots/normIa_resid_v_mass.png')
                 # title='Hubble Residual v. Host Stellar Mass of Normal SNe Ia from CSP [SNooPy]',)

    # Update Redshift Plots
    resid_v_z(path='output/combiend__snpy_params_cut.txt',
              save_loc='saved/readme_plots/csp-atlas-ztf_snpy_resid_v_z.png')
              # title='Hubble Residual v. CMB Redshift of CSP-ATLAS-ZTF 91bg-like SNe Ia [SNooPy]')
    resid_v_z(path='output/combiend__salt_params_cut.txt',
              save_loc='saved/readme_plots/csp-atlas-ztf_salt_resid_v_z.png')
              # title = 'Hubble Residual v. CMB Redshift of CSP-ATLAS-ZTF 91bg-like SNe Ia [SALT3]',
    resid_v_z(path='output/merged_params_cut.txt',
              save_loc='saved/readme_plots/merged_resid_v_z.png')
              # title = 'Hubble Residual v. CMB Redshift of CSP-ATLAS-ZTF 91bg-like SNe Ia [SALT3-SNooPy]',

    # Update Histograms
    param_hist('output/combiend__snpy_params_cut.txt',
               'output/norm_snpy_params_cut.txt',
               algo='snpy', type='median', st_width=0.04, c_width=0.04,
               save_loc='saved/readme_plots/snpy_params_hicat_v_dr3.png')
    param_hist('output/combiend__salt_params_cut.txt',
               'output/norm_salt_params_cut.txt',
               algo='salt', type='median', st_width=0.3, c_width=0.08,
               save_loc='saved/readme_plots/salt_params_hicat_v_dr3.png')
    return

if __name__ == '__main__':
    start = systime.time()  # Runtime tracker

    # format_dr3()

    refresh_all()
    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')
