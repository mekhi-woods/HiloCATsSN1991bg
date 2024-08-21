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

class sn91bg_snpy():
    def __init__(self, objname=None, coords=(0.00, 0.00), z=0.00, origin=None):
        self.objname = objname
        self.coords = coords
        self.z = z
        self.origin = origin

        self.period = (999999.9, 999999.9)
        self.mu = {'value': 0.00, 'err': 0.00}
        self.st = {'value': 0.00, 'err': 0.00}
        self.Tmax = {'value': 0.00, 'err': 0.00}
        self.EBVhost = {'value': 0.00, 'err': 0.00}
        self.hostMass = {'value': 0.00, 'err': 0.00}

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
                    '---------------------------------------------------------------------------------------------\n' +
                    '\tStretch (st) = ' + str(self.st['value']) + ' +/- ' + str(self.st['err']) + '\n' +
                    '\tDistance Mod. (mu) = ' + str(self.mu['value']) + ' +/- ' + str(self.mu['err']) + '\n' +
                    '\tTmax = ' + str(self.Tmax['value']) + ' +/- ' + str(self.Tmax['err']) + '\n' +
                    '\tEBVHost = ' + str(self.EBVhost['value']) + ' +/- ' + str(self.EBVhost['err']) + '\n' +
                    '\tHost Mass = ' + str(self.hostMass['value']) + ' +/- ' + str(self.hostMass['err']) + '\n' +
                    '---------------------------------------------------------------------------------------------')
        print(prnt_str)
    def save_class(self, save_loc):
        print('[+++] '+self.objname+' -- Saving class to '+save_loc+'classes/'+self.objname+'_class.txt')
        with open(save_loc+'classes/'+self.objname+'_class.txt', 'w') as f:
            f.write(self.origin + ' ' + self.objname + ' ' + str(self.coords[0]) + ' ' + str(self.coords[1]) + ' ' + str(self.z) + '\n')
            f.write('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
            f.write('st, ' + str(self.st['value']) + ', ' + str(self.st['err']) + '\n')
            f.write('mu, ' + str(self.mu['value']) + ', ' + str(self.mu['err']) + '\n')
            f.write('Tmax, ' + str(self.Tmax['value']) + ', ' + str(self.Tmax['err']) + '\n')
            f.write('EBVHost, ' + str(self.EBVhost['value']) + ', ' + str(self.EBVhost['err']) + '\n')
            f.write('HostMass, ' + str(self.hostMass['value']) + ', ' + str(self.hostMass['err']) + '\n')
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

            for cat in [self.st, self.mu, self.Tmax, self.EBVhost, self.hostMass]:
                line = f.readline().split(', ')[1:]
                cat['value'], cat['err'] = float(line[0]), float(line[1])

            skip = f.readline()

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
        print('Best filters:', list(n_s.data.keys()))

        for i in range(5):
            try:
                n_s.fit(bands=None, dokcorr=True, k_stretch=False, reset_kcorrs=True, **{'mangle': 1, 'calibration': 0})
                n_s.save(save_path)

                self.mu = {'value': n_s.parameters['DM'], 'err': n_s.errors['DM']}
                self.st = {'value': n_s.parameters['st'], 'err': n_s.errors['st']}
                self.Tmax = {'value': n_s.parameters['Tmax'], 'err': n_s.errors['Tmax']}
                self.EBVhost = {'value': n_s.parameters['EBVhost'], 'err': n_s.errors['EBVhost']}

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
                    break
                else:
                    print(error)

        # Restore print statements
        sys.stdout = sys.__stdout__

        return
class sn91bg_salt():
    def __init__(self, objname=None, coords=(0.00, 0.00), z=0.00, origin=None):
        self.objname = objname
        self.coords = coords
        self.z = z
        self.origin = origin

        self.period = (999999.9, 999999.9)
        self.mu = {'value': 0.00, 'err': 0.00}
        self.t0 = {'value': 0.00, 'err': 0.00}
        self.x0 = {'value': 0.00, 'err': 0.00}
        self.x1 = {'value': 0.00, 'err': 0.00}
        self.c = {'value': 0.00, 'err': 0.00}
        self.hostMass = {'value': 0.00, 'err': 0.00}

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
                    '---------------------------------------------------------------------------------------------\n' +
                    '\tDistance Mod. (mu) = ' + str(self.mu['value']) + ' +/- ' + str(self.mu['err']) + '\n' +
                    '\tT-max (t0) = ' + str(self.t0['value']) + ' +/- ' + str(self.t0['err']) + '\n' +
                    '\tAmplitude (x0) = ' + str(self.x0['value']) + ' +/- ' + str(self.x0['err']) + '\n' +
                    '\tTime-dependent Variation (x1) = ' + str(self.x1['value']) + ' +/- ' + str(self.x1['err']) + '\n' +
                    '\tColor (c) = ' + str(self.c['value']) + ' +/- ' + str(self.c['err']) + '\n' +
                    '\tHost Mass = ' + str(self.hostMass['value']) + ' +/- ' + str(self.hostMass['err']) + '\n' +
                    '---------------------------------------------------------------------------------------------')
        print(prnt_str)
    def save_class(self, save_loc):
        print('[+++] '+self.objname+' -- Saving class to '+save_loc+'classes/'+self.objname+'_class.txt')
        with open(save_loc+'classes/'+self.objname+'_class.txt', 'w') as f:
            f.write(self.origin + ' ' + self.objname + ' ' + str(self.coords[0]) + ' ' + str(self.coords[1]) + ' ' + str(self.z) + '\n')
            f.write('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
            f.write('mu, ' + str(self.mu['value']) + ', ' + str(self.mu['err']) + '\n')
            f.write('t0, ' + str(self.t0['value']) + ', ' + str(self.t0['err']) + '\n')
            f.write('x0, ' + str(self.x0['value']) + ', ' + str(self.x0['err']) + '\n')
            f.write('x1, ' + str(self.x1['value']) + ', ' + str(self.x1['err']) + '\n')
            f.write('c, ' + str(self.c['value']) + ', ' + str(self.c['err']) + '\n')
            f.write('HostMass, ' + str(self.hostMass['value']) + ', ' + str(self.hostMass['err']) + '\n')
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

            for cat in [self.mu, self.t0, self.x0, self.x1, self.c, self.hostMass]:
                line = f.readline().split(', ')[1:]
                cat['value'], cat['err'] = float(line[0]), float(line[1])

            skip = f.readline()

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

            self.t0['value'], self.t0['err'] = result.parameters[1], result.errors['t0']
            self.x0['value'], self.x0['err'] = result.parameters[2], result.errors['x0']
            self.x1['value'], self.x1['err'] = result.parameters[3], result.errors['x1']
            self.c['value'], self.c['err'] = result.parameters[4], result.errors['c']


            # Calculate
            pho_mB = -2.5 * np.log10(self.x0['value']) + mB_const
            pho_mB_err = np.abs((-2.5 * self.x0['err']) / (self.x0['err'] * np.log(10)))

            mu = pho_mB + (alpha * self.x1['value']) - (beta * self.c['value']) - M0
            mu_err = np.sqrt(pho_mB_err ** 2 + (np.abs(alpha) * self.x1['err']) ** 2 + (
                        np.abs(alpha) * self.c['err']) ** 2)

            self.mu['value'], self.mu['err'] = mu, mu_err

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
            self.mu['value'] = -1.0

        return
def snpy_class_creation(data_set, path, save=True, results=False):
    if data_set == 'CSP':
        with open(path, 'r') as f:
            objname, z, ra, dec = f.readline().split(' ')
        objname, ra, dec, z = objname[2:], float(ra), float(dec[:-1]), float(z)
        tempSN = sn91bg_snpy(objname, (ra, dec), z, 'CSP')

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
    elif data_set == 'ATLAS':
        data = np.genfromtxt(path, delimiter=',', dtype=str, skip_header=1)
        if len(data) == 0:
            print('[!!!] File [' + path + '] empty!')
            return
        ra, dec = np.average(data[:, 1].astype(float)), np.average(data[:, 2].astype(float))
        objname, z = gen.TNS_objname_z(ra, dec)
        z = np.nan if z == 'None' else float(z)
        tempSN = sn91bg_snpy(objname, (ra, dec), z, 'ATLAS')
        tempSN.set_data(zp=data[:, 7].astype(float), filters=data[:, 6], time=data[:, 8],
                        flux=data[:, 16], dflux=data[:, 17], mag=data[:, 3], dmag=data[:, 4])
    elif data_set == 'ZTF':
        data = np.genfromtxt(path, delimiter=None, dtype=str, skip_header=56)
        if len(data) == 0:
            print('[!!!] File [' + path + '] empty!')
            return

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
            return

        tempSN = sn91bg_snpy(objname, (ra, dec), z, 'ZTF')
        tempSN.set_data(zp=zp, filters=filters, time=time,
                        flux=flux, dflux=dflux, mag=mag, dmag=dmag)
    else:
        raise ValueError("Data set '" + data_set + "' not recognized")

    tempSN.clean_data()
    tempSN.write_snpy_ASCII(save_loc=CONSTANTS[data_set.lower() + '_saved_loc'] + 'ascii/')
    tempSN.snpy_fit(save_loc=CONSTANTS[data_set.lower() + '_saved_loc'])
    # tempSN.get_host_mass()

    if save:
        tempSN.save_class(save_loc=CONSTANTS[data_set.lower() + '_saved_loc'])
    if results:
        tempSN.print_info()
    return tempSN
def snpy_batch(data_set, save=True, results=False):
    SNe = []
    if data_set == 'COMBINED':
        for set in ['CSP', 'ATLAS', 'ZTF']:
            files = glob.glob('../data/' + set + '/*.txt')
            for path in files:
                print('[', files.index(path) + 1, '/', len(files), ']')
                SNe.append(snpy_class_creation(set, path, save=True, results=True))
    else:
        files = glob.glob('../data/' + data_set + '/*.txt')
        for path in files:
            print('[', files.index(path) + 1, '/', len(files), ']')
            SNe.append(snpy_class_creation(data_set, path, save=True, results=True))
    return SNe
def salt_class_creation(data_set, path, save=True, results=False):
    if data_set == 'CSP':
        with open(path, 'r') as f:
            objname, z, ra, dec = f.readline().split(' ')
        objname, ra, dec, z = objname[2:], float(ra), float(dec[:-1]), float(z)
        tempSN = sn91bg_salt(objname, (ra, dec), z, 'CSP')

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
    elif data_set == 'ATLAS':
        data = np.genfromtxt(path, delimiter=',', dtype=str, skip_header=1)
        if len(data) == 0:
            print('[!!!] File [' + path + '] empty!')
            return
        ra, dec = np.average(data[:, 1].astype(float)), np.average(data[:, 2].astype(float))
        objname, z = gen.TNS_objname_z(ra, dec)
        z = np.nan if z == 'None' else float(z)
        tempSN = sn91bg_salt(objname, (ra, dec), z, 'ATLAS')
        tempSN.set_data(zp=data[:, 7].astype(float), filters=data[:, 6], time=data[:, 8],
                        flux=data[:, 16], dflux=data[:, 17], mag=data[:, 3], dmag=data[:, 4])
    elif data_set == 'ZTF':
        data = np.genfromtxt(path, delimiter=None, dtype=str, skip_header=56)
        if len(data) == 0:
            print('[!!!] File [' + path + '] empty!')
            return

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
            return

        tempSN = sn91bg_salt(objname, (ra, dec), z, 'ZTF')
        tempSN.set_data(zp=zp, filters=filters, time=time,
                        flux=flux, dflux=dflux, mag=mag, dmag=dmag)
    else:
        raise ValueError("Data set '" + data_set + "' not recognized")

    tempSN.clean_data()
    tempSN.salt_fit(save_loc=CONSTANTS['salt_'+data_set.lower()+'_loc'])
    # tempSN.get_host_mass()

    if save:
        tempSN.save_class(save_loc=CONSTANTS[data_set.lower() + '_saved_loc'])
    if results:
        tempSN.print_info()
    return tempSN
def salt_batch(data_set, save=True, results=False, readout=False, pause_error=False):
    SNe = []
    if data_set == 'COMBINED':
        for set in ['CSP', 'ATLAS', 'ZTF']:
            files = glob.glob('../data/' + set + '/*.txt')
            for path in files:
                print('[', files.index(path) + 1, '/', len(files), ']')
                SNe.append(salt_class_creation(set, path, save=True, results=True))
    else:
        files = glob.glob('../data/' + data_set + '/*.txt')
        for path in files:
            print('[', files.index(path) + 1, '/', len(files), ']')
            tempSN = salt_class_creation(data_set, path, save=True, results=True)
            if tempSN.mu['value'] <= 0:
                print('[!!!] SALT fit failed, excluding...')
                if pause_error :
                    input('[!!!] An error has occured, press any key to continue...')
                continue
            else:
                SNe.append(tempSN)
    if readout:
        for SN in SNe:
            print(SN)
            print('\tmu:', SN.mu['value'], '+/-', SN.mu['err'], 't0:', SN.t0['value'], '+/-', SN.t0['err'])
            print('\tx0:', SN.x0['value'], '+/-', SN.x0['err'], 'x1:', SN.x1['value'], '+/-', SN.x1['err'])
            print('\tc:', SN.c['value'], '+/-', SN.c['err'])
            print('----------------------------------------------------------------')
    return SNe
def batch_load_from_file(data_set):
    SNe = []
    if data_set == 'COMBINED':
        for set in ['CSP', 'ATLAS', 'ZTF']:
            files = glob.glob('../saved/snpy/' + set.lower() + 'classes/*_class.txt')
            for path in files:
                print('[', files.index(path) + 1, '/', len(files), ']')
                tempSN = sn91bg_snpy()
                tempSN.load_from_file(path)
                SNe.append(tempSN)
    else:
        files = glob.glob('../saved/snpy/' + data_set.lower() + 'classes/*_class.txt')
        for path in files:
            print('[', files.index(path) + 1, '/', len(files), ']')
            tempSN = sn91bg_snpy()
            tempSN.load_from_file(path)
            SNe.append(tempSN)
    return SNe

if __name__ == '__main__':
    # Load and initiallize in batch
    # CSPSNe = batch_class_creation('CSP', save=True, results=False)

    # Load and initiallize indivisual
    # SN2006bd = class_creation('CSP', '../data/CSP/SN2006bd_snpy.txt', save=True, results=False)
    # SN2006bd.print_info()

    # Load from file indivisual
    # SN2006bd = sn91bg()
    # SN2006bd.load_from_file('../saved/snpy/csp/classes/2006bd_class.txt')
    # SN2006bd.print_info()

    # allSNe = batch_class_creation('COMBINED', save=True, results=True)

    # +============================================================+ #

    # SN2006bd = salt_class_creation('CSP', '../data/CSP/SN2006bd_snpy.txt', save=True, results=False)

    salt_batch('ATLAS', save=True, results=True, readout=True, pause_error=True)
