# def null():
#     import warnings
#     with warnings.catch_warnings():
#         import scripts.general as gen
#         from astropy.cosmology import FlatLambdaCDM
#         import matplotlib.pyplot as plt
#         import numpy as np
#         import sys
#         import os
#
#     CURRENT_COSMO = FlatLambdaCDM(70, 0.3)  # Hubble Constant, Omega-Matter
#     CONSTANTS = gen.get_constants()
#
#     def norm_Ia_objs_pull(loc='../txts/normIaparams.txt', quiet=False):
#         # Check quiet
#         if quiet:
#             sys.stdout = open(os.devnull, 'w')
#
#         normIaObjs = {}
#         with open(loc) as f:
#             hdr = []
#             for n in f.readline().split(' '):
#                 if n != '':
#                     hdr.append(n)
#             hdr[-1] = hdr[-1][:-1]
#             print(hdr)
#
#             data = f.readlines()
#             i, t = 0, len(data[:-1])
#             for line in data[:-1]:
#                 i += 1
#                 print(
#                     '---------------------------------------------------------------------------------------------------')
#                 print('[ ' + str(i) + ' / ' + str(t) + ' ]')
#                 n_line = line.split(' ')
#                 n_line[-1] = n_line[-1][:-1]  # Remove '\n' from last element
#                 n_line[0] = n_line[0][:-4]  # Remove the .txt from ID
#                 ra, dec = float(n_line[1]), float(n_line[2])  # Save RA/DEC
#
#                 objname, z = gen.TNS_objname_z(ra, dec)  # Save Objname
#
#                 normIaObjs.update({objname: {
#                     'ra': ra, 'dec': dec, 'z': z, 'mu': n_line[4], 'mu_err': 0.00,
#                     't0': n_line[5], 'x0': n_line[6], 'x1': n_line[7], 'c': n_line[8],
#                     't0err': n_line[9], 'x0err': n_line[10], 'x1err': n_line[11], 'cerr': n_line[12]
#                     # 'original_name': n_line[0],
#                 }})
#
#         # Restore print statements
#         sys.stdout = sys.__stdout__
#
#         return normIaObjs
#
#     def sort_resid_z(objs, x1lim=1):
#         # objnames, all_resid, all_resid_err, all_z = [], [], [], []
#         # for obj in objs:
#         #     # Normalize Data
#         #     if str(objs[obj]['mu']) == 'nan' or str(objs[obj]['z']) == 'nan':
#         #         continue
#         #     elif str(objs[obj]['mu']) == 'None' or str(objs[obj]['z']) == 'None':
#         #         continue
#         #     # Last adjustments
#         #     elif float(objs[obj]['x1']) < -x1lim or float(objs[obj]['x1']) > x1lim:
#         #         continue
#         #     else:
#         #         z = float(objs[obj]['z'])
#         #         resid = float(objs[obj]['mu']) - CURRENT_COSMO.distmod(float(objs[obj]['z'])).value
#         #         resid_err = float(objs[obj]['mu_err'])
#         #
#         #         all_z.append(z)
#         #         all_resid.append(resid)
#         #         all_resid_err.append(resid_err)
#         #         objnames.append(obj)
#         # return objnames, all_resid, all_resid_err, all_z
#
#         return
#
#     def plot_z_v_mu(objnames, all_resid, all_resid_err, all_z, title, x1lim='', labels=False):
#         fig, axs = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'width_ratios': [10, 1]}, constrained_layout=True)
#
#         # Plotting
#         axs[0].errorbar(all_z, all_resid, yerr=all_resid_err, fmt='o')
#         axs[1].hist(all_resid, bins=40, orientation="horizontal")
#         if labels:
#             for i in range(len(objnames)):
#                 axs[0].text(all_z[i], all_resid[i], objnames[i], fontsize=8)
#
#         # Formatting
#         ylimiter = np.max(np.abs(all_resid)) + 0.5
#         axs[0].set_ylim(-ylimiter, ylimiter);
#         axs[1].set_ylim(-ylimiter, ylimiter)
#         fig.suptitle(title +  # Figure Title
#                      '-' + str(x1lim) + ' < x1 < ' + str(x1lim) + ' | Median: ' + str(
#             round(np.median(all_resid), 2)) + ' | '
#                                               'Scatter: ' + str(round(np.std(all_resid), 2)) + ' | # of pts: ' + str(
#             len(all_resid)), size='medium')
#         axs[0].set(xlabel='Redshift', ylabel='Hubble Residuals')  # Sub-plot Labels
#         axs[1].get_yaxis().set_visible(False)  # Turn off y-axis labels
#         plt.show()
#         return
#
#     def plot_mass_v_mu(objnames, all_resid, all_z, title, x1lim='', labels=False):
#         fig, axs = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'width_ratios': [10, 1]}, constrained_layout=True)
#
#         # Plotting
#         axs[0].scatter(all_z, all_resid)
#         axs[1].hist(all_resid, bins=40, orientation="horizontal")
#         if labels:
#             for i in range(len(objnames)):
#                 axs[0].text(all_z[i], all_resid[i], objnames[i], fontsize=8)
#
#         # Formatting
#         ylimiter = np.max(np.abs(all_resid)) + 0.5
#         axs[0].set_ylim(-ylimiter, ylimiter);
#         axs[1].set_ylim(-ylimiter, ylimiter)
#         fig.suptitle(title +  # Figure Title
#                      '-' + str(x1lim) + ' < x1 < ' + str(x1lim) + ' | Median: ' + str(
#             round(np.median(all_resid), 2)) + ' | '
#                                               'Scatter: ' + str(round(np.std(all_resid), 2)) + ' | # of pts: ' + str(
#             len(all_resid)), size='medium')
#         axs[0].set(xlabel='Redshift', ylabel='Hubble Residuals')  # Sub-plot Labels
#         axs[1].get_yaxis().set_visible(False)  # Turn off y-axis labels
#         plt.show()
#         return
#
#     if __name__ == '__main__':
#         # params = {'x1lim': 0.5}
#         #
#         # # Plot Normal SNe Ia
#         # # txtNormPath = '../default/sneNormObjs.txt'
#         #
#         # # sneNormObjs = norm_Ia_objs_pull(quiet=True)
#         # # objnames, all_resid, all_resid_err, all_z = sort_resid_z(sneNormObjs, **params)
#         # # gen.dict_handler(choice='pack', data_dict=sneNormObjs, path=txtNormPath)
#         #
#         # # gen.host_mass(txtNormPath, save_loc='../default/', keep_data=False, update_saved=True, use_mass_key=True)
#         #
#         # # plot_z_v_mu(objnames, all_resid, all_resid_err, all_z, title='Hubble Residuals vs. Redshift of Normal SNe Ia\n', labels=True, **params)
#         #
#         #
#         # # Plot 91bg-like SNe Ia
#         # # txt91bgPath = '../default/sne91bgObjs.txt'
#         # #
#         # # sn91bgIaObjs = gen.dict_handler(choice='unpack', path='../saved/salt/atlas/salt_atlas_saved.txt')
#         # # print(len(sn91bgIaObjs))
#         # # objnames, all_resid, all_resid_err, all_z = sort_resid_z(sn91bgIaObjs, **params)
#         # #
#         # # gen.dict_handler(choice='pack', data_dict=sn91bgIaObjs, path=txt91bgPath)
#         # gen.host_mass(txt91bgPath, save_loc='../default/', keep_data=False, update_saved=True, use_mass_key=True)
#         #
#         # plot_z_v_mu(objnames, all_resid, all_resid_err, all_z, title='Hubble Residuals vs. Redshift of 91bg-like SNe Ia\n', labels=True, **params)
#
#         # gen.mass_step_calc(txtNormPath)
#
#         objs = gen.dict_handler(choice='unpack', path=CONSTANTS['atlas_saved_loc'] + 'atlas_saved.txt')
#
#     return

import warnings
warnings.simplefilter("ignore", UserWarning)
import sys
import snpy
import glob
import sncosmo
import numpy as np
import time as systime
import matplotlib.pyplot as plt
from astropy.table import Table

import scripts.general as gen

CONSTANTS = gen.get_constants()

class sn91bg():
    def __init__(self, objname=None, coords=(0.00, 0.00), z=0.00, origin=None):
        self.objname = objname
        self.coords = coords
        self.z = z
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
        return self.origin + ' | ' + self.objname + ' | ' + str(self.coords) + ' | z = ' + str(self.z)
    def print_info(self):
        prnt_str = ('---------------------------------------------------------------------------------------------\n' +
                    self.origin + ' | ' + self.objname + ' | ' + str(self.coords) + ' | z = ' + str(self.z) + '\n' +
                    '---------------------------------------------------------------------------------------------\n' +
                    '\tFilters = ' + str(np.unique(self.filters)) + '\n' +
                    '\tFlux Range = (' + str(np.min(self.flux)) + ' +/- ' + str(self.dflux[np.where(self.flux == np.min(self.flux))[0][0]]) +
                    ') -- (' + str(np.max(self.flux)) + ' +/- ' + str(self.dflux[np.where(self.flux == np.max(self.flux))[0][0]]) + ')\n' +
                    '\tMagnitude Range = (' + str(np.min(self.mag)) + ' +/- ' + str(self.dmag[np.where(self.mag == np.min(self.mag))[0][0]]) +
                    ') -- (' + str(np.max(self.mag)) + ' +/- ' + str(self.dmag[np.where(self.mag == np.max(self.mag))[0][0]]) + ')\n' +
                    '\tMJDs-MJDe = ' + str(self.period[0]) + ' -- ' + str(self.period[1]) + '\n' +
                    '---------------------------------------------------------------------------------------------\n')
        for p in self.params:
            prnt_str += ('\t'+p+' = ' + str(self.params[p]['value']) + ' +/- ' + str(self.params[p]['err']) + '\n')
        prnt_str += '---------------------------------------------------------------------------------------------'
        print(prnt_str)
    def save_class(self, save_loc):
        print('[+++] '+self.objname+' -- Saving class to '+save_loc+'classes/'+self.objname+'_class.txt')
        with open(save_loc+'classes/'+self.objname+'_class.txt', 'w') as f:
            f.write(self.origin + ' ' + self.objname + ' ' + str(self.coords[0]) + ' ' + str(self.coords[1]) + ' ' + str(self.z) + '\n')
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
            self.origin, self.objname, self.coords, self.z = hdr[0], hdr[1], (float(hdr[2]), float(hdr[3])), float(hdr[4][:-1])

            skip = f.readline()

            # for cat in [self.mu, self.t0, self.x0, self.x1, self.c, self.hostMass]:
            #     line = f.readline().split(', ')[1:]
            #     cat['value'], cat['err'] = float(line[0]), float(line[1])

            line = f.readline()
            while '+++' not in line:
                line = line.split(', ')
                self.params.update({line[0]: {'value': float(line[1]), 'err': float(line[2])}})
                line = f.readline()

            self.zp = np.array(f.readline().split(', ')[1:]).astype(float)
            self.filters = np.array(f.readline().split(', ')[1:]).astype(str)
            self.time = np.array(f.readline().split(', ')[1:]).astype(float)
            self.flux = np.array(f.readline().split(', ')[1:]).astype(float)
            self.dflux = np.array(f.readline().split(', ')[1:]).astype(float)
            self.mag = np.array(f.readline().split(', ')[1:]).astype(float)
            self.dmag = np.array(f.readline().split(', ')[1:]).astype(float)
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
    def clean_data(self, mag_unc_max=0, flux_unc_max=0):
        print('[+++] '+self.objname+' -- Cleaning data...')
        new_zp, new_filters, new_time, new_mag, new_dmag, new_flux, new_dflux = (
            np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]))
        for n in range(len(self.zp)):
            n_zp, n_filters, n_time  = self.zp[n], self.filters[n], self.time[n]
            n_mag, n_dmag, n_flux, n_dflux = self.mag[n], self.dmag[n], self.flux[n], self.dflux[n]

            n_mag = str(n_mag).replace('>', '')
            if n_mag == 'None' or n_dmag == 'None' or n_flux == 'None' or n_dflux == 'None':
                continue
            if float(n_dflux) == 0 or float(n_dmag) == 0:
                continue
            if float(n_mag) <= 0 or float(n_flux) <= 0:
                continue
            if (mag_unc_max != 0) and (float(n_dmag) > mag_unc_max):
                continue
            if (flux_unc_max != 0) and (float(n_dflux) > flux_unc_max):
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
    def write_snpy_ASCII(self, save_loc='../default/'):
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
            n_s.k_version = '91bg'
            n_s.choose_model('EBV_model2', stype='st')
            n_s.set_restbands()  # Auto pick appropriate rest-bands
        except:
            print('[!!!] Failed to load ASCII file')
            return

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
            salt_filters = np.array([])
            for filter in self.filters:
                salt_filters = np.append(salt_filters, filter_dict[filter])

            data = Table([self.time, salt_filters, self.flux, self.dflux,
                          self.zp, np.full(len(self.time), 'ab')],
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

            print('Pausing for 1 seconds...')
            systime.sleep(1)
        except Exception as error:
            print(error)
            self.params.update({'mu': {'value': -1.0, 'err': -1.0}})

        return
    def get_host_mass(self):
        return

def class_creation(data_set, path):
    if data_set == 'CSP':
        with open(path, 'r') as f:
            objname, z, ra, dec = f.readline().split(' ')
        objname, ra, dec, z = objname[2:], float(ra), float(dec[:-1]), float(z)
        tempSN = sn91bg(objname, (ra, dec), z, 'CSP')

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
                    n_flux = 10 ** ((n_mag - zp) / 2.5)
                    n_dflux = 10 ** (n_dmag / 2.5)

                    tempSN.zp = np.append(tempSN.zp, zp)
                    tempSN.filters = np.append(tempSN.filters, filter)
                    tempSN.time = np.append(tempSN.time, n_time)
                    tempSN.flux = np.append(tempSN.flux, n_flux)
                    tempSN.dflux = np.append(tempSN.dflux, n_dflux)
                    tempSN.mag = np.append(tempSN.mag, n_mag)
                    tempSN.dmag = np.append(tempSN.dmag, n_dmag)

        tempSN.clean_data()
    elif data_set == 'ATLAS':
        data = np.genfromtxt(path, delimiter=',', dtype=str, skip_header=1)
        if len(data) == 0:
            print('[!!!] File [' + path + '] empty!')
            return None
        ra, dec = np.average(data[:, 1].astype(float)), np.average(data[:, 2].astype(float))
        objname, z = gen.TNS_objname_z(ra, dec)
        z = np.nan if z == 'None' else float(z)
        tempSN = sn91bg(objname, (ra, dec), z, 'ATLAS')
        tempSN.set_data(zp=data[:, 7].astype(float), filters=data[:, 6], time=data[:, 8],
                        flux=data[:, 16], dflux=data[:, 17], mag=data[:, 3], dmag=data[:, 4])
        tempSN.clean_data()
    elif data_set == 'ZTF':
        data = np.genfromtxt(path, delimiter=None, dtype=str, skip_header=56)
        if len(data) == 0:
            print('[!!!] File [' + path + '] empty!')
            return None
        with (open(path, 'r') as f):
            hdr = f.readlines()
            ra, dec = float(hdr[3].split(' ')[-2]), float(hdr[4].split(' ')[-2])
            objname, z = gen.TNS_objname_z(ra, dec)
            z = np.nan if z == 'None' else float(z)

            # Get magnitudes m = -2.5log(F) + zp
            time, flux, dflux, zp, filters = data[:, 22], data[:, 24], data[:, 25], data[:, 20], data[:, 4]
            valid_ints = np.unique(np.hstack((np.where(flux != 'null')[0], np.where(dflux != 'null')[0])))
            time, zp, filters = time[valid_ints].astype(float), zp[valid_ints].astype(float), filters[valid_ints]
            flux, dflux = flux[valid_ints].astype(float), dflux[valid_ints].astype(float)
            mag, dmag = ((-2.5 * np.log10(flux)) + zp), (2.5 * np.log10(dflux))

        if len(time) == 0:
            print('[!!!] File has no valid values!')
            return None

        tempSN = sn91bg(objname, (ra, dec), z, 'ZTF')
        tempSN.set_data(zp=zp, filters=filters, time=time,
                        flux=flux, dflux=dflux, mag=mag, dmag=dmag)
        tempSN.clean_data()
    else:
        raise ValueError("Data set '" + data_set + "' not recognized")
    return tempSN
def indvisual_fit(data_set, path, algo='snpy'):
    if algo == 'snpy':
        tempSN = class_creation(data_set, path)
        if tempSN is None:
            return None
        tempSN.write_snpy_ASCII(save_loc=CONSTANTS[data_set.lower() + '_saved_loc'] + 'ascii/')
        tempSN.snpy_fit(save_loc=CONSTANTS[data_set.lower() + '_saved_loc'])
        if tempSN.params['mu']['value'] <= 0.00:
            return None
        tempSN.get_host_mass()
        tempSN.save_class(CONSTANTS[data_set.lower() + '_saved_loc'])
    elif algo == 'salt':
        tempSN = class_creation(data_set, path)
        if tempSN is None:
            return None
        tempSN.salt_fit(save_loc=CONSTANTS['salt_' + data_set.lower() + '_loc'])
        if tempSN.params['mu']['value'] <= 0.00:
            return None
        tempSN.get_host_mass()
        tempSN.save_class(CONSTANTS['salt_' + data_set.lower() + '_loc'])
    else:
        raise ValueError("[!!!] Invalid algorithm selected ['snpy'/'salt']")
    return tempSN
def batch_fit(data_set, algo='snpy'):
    SNe = []
    num_int = 0
    if data_set == 'COMBINED':
        for d_set in ['CSP', 'ATLAS', 'ZTF']:
            files = glob.glob('../data/' + d_set + '/*.txt')
            num_int += len(files)
            for path in files:
                print('[', files.index(path) + 1, '/', len(files), ']')
                print('-----------------------------------------------------------------------------------------------')

                tempSN = indvisual_fit(d_set, path, algo)
                if tempSN is not None:
                    SNe.append(tempSN)
    else:
        files = glob.glob('../data/' + data_set + '/*.txt')
        num_int += len(files)
        for path in files:
            print('[', files.index(path) + 1, '/', len(files), ']')
            print('-----------------------------------------------------------------------------------------------')
            tempSN = indvisual_fit(data_set, path, algo)
            if tempSN is not None:
                SNe.append(tempSN)
    print('Sucessfully fit [', len(SNe), '/', num_int, ']!')

    # Save params to file
    if algo == 'snpy':
        save_loc = CONSTANTS[data_set.lower() + '_saved_loc']
    if algo == 'salt':
        save_loc = CONSTANTS['salt_'+data_set.lower()+'_loc']
    with open(save_loc + data_set.lower() + '_params.txt', 'w') as f:
        hdr = 'objname, ra, dec, z, MJDs, MJDe, origin'
        for params in SNe[0].params:
            hdr += ', ' + str(params) + ', ' + str(params) + '_err'
        f.write(hdr + '\n')
        for SN in SNe:
            line = (SN.objname + ', ' + str(SN.coords[0]) + ', ' + str(SN.coords[1]) + ', ' + str(SN.z) + ', ' +
                    str(SN.period[0]) + ', ' + str(SN.period[1]) + ', ' + str(SN.origin))
            for param in SN.params:
                line += ', ' + str(SN.params[param]['value']) + ', ' + str(SN.params[param]['err'])
            f.write(line + '\n')
    return SNe
def indivisual_load():
    return
def batch_load():
    return
if __name__ == '__main__':
    start = systime.time() # Runtime tracker

    # SN = indvisual_fit('CSP', '../data/CSP/SN2006bd_snpy.txt', algo='snpy')
    SNe = batch_fit('CSP', 'snpy')

    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')
