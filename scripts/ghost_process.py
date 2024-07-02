from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from astropy import cosmology as cosmo
from astro_ghost.ghostHelperFunctions import getTransientHosts
import shutil
import numpy as np

import matplotlib.pyplot as plt

import scripts.general as gen

COSMO_MODEL = cosmo.FlatLambdaCDM(H0=70, Om0=0.3)

GHOST_DATA = '../data/GHOST/'
TEST_ROOT = '../tests/'

ATLAS_SAVE_TXT = '../snpy/atlas/atlas_saved.txt'
BURNS_SAVE_TXT = '../snpy/burns/burns_saved.txt'

def ghost_host_galaxy(dict_path, save_loc=TEST_ROOT, keep_data=True, update_saved=False):
    cosmo = FlatLambdaCDM(70, 0.3) # Hubble Constant, Omega-Matter
    data = gen.dict_unpacker(dict_path)
    all_z, all_logstellarmass = [], []
    print('[+++] Finding host galaxy mass using GHOST...')

    i = 0
    for obj in data:
        ra, dec, z = float(data[obj]['ra']), float(data[obj]['dec']), float(data[obj]['z'])

        print('\n[', list(data).index(obj)+1, '/', len(data), ']', obj, '|', data[obj]['ra'], data[obj]['dec'], data[obj]['z'])
        print('---------------------------------------------------------------------------------------------------------')

        transient_position = SkyCoord(ra, dec, unit=u.deg)
        try:
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                host_data = getTransientHosts(transientCoord=[transient_position], transientName=[obj], verbose=False, starcut="gentle", savepath=save_loc+'ghost_stuff/', GHOSTpath=GHOST_DATA)

            gMag, iMag, iAbsMag = host_data['gKronMag'].loc[0], host_data['iKronMag'].loc[0], (host_data['iKronMag'].loc[0] - cosmo.distmod(z).value)
            gMagErr, iMagErr, iAbsMagErr = host_data['gKronMagErr'].loc[0], host_data['iKronMagErr'].loc[0], host_data['iKronMagErr'].loc[0]
            giMagErr = np.sqrt((gMagErr**2) + (iMagErr**2))

            logstellarmass = (1.15 + (0.7*(gMag - iMag)) - (0.4*(iAbsMag))) # Taylor et. al. 2011 -- eq. 8
            logstellarmasserr = np.sqrt(((0.7**2)*(giMagErr**2)) + ((0.4**2)*(iAbsMagErr**2)))

        except Exception as error:
            print(obj+':', error)
            logstellarmass, logstellarmasserr = 0.00, 0.00

        if logstellarmass != np.nan and logstellarmass > 0:
            print('Success!', obj, 'host galaxy has a mass of:', logstellarmass, '+/-', logstellarmasserr, 'logM_* / [Msun]')
            all_z.append(z)
            all_logstellarmass.append(logstellarmass)
            if update_saved:
                data[obj].update({'logstellarmass': logstellarmass, 'logstellarmasserr': logstellarmasserr})
        else:
            print('Failed to find host galaxy!')
            if update_saved:
                data[obj].update({'logstellarmass': 0.00, 'logstellarmasserr': 0.00})

    print('\nSuccessfully found mass of', len(all_z), '/', len(data), 'host galaxies!')
    if not keep_data:
        print('Removing GHOST data...')
        shutil.rmtree(save_loc+'ghost_stuff/') # Clear messy data
    if update_saved:
        print('Saving data to'+dict_path+'...')
        gen.dict_packer(data, dict_path)

    return all_z, all_logstellarmass
def ghost_plotting(choice, plot_size = (18, 6), plot_ratio = [10, 1], hist_bins = [50, 50, 10], labels = True, raw = False):
    if 'atlas_all' in choice:
        choice = ['altas_muvmass', 'altas_muvz', 'altas_muvmu', 'altas_residualsvz', 'altas_residualsvmass']
    if 'burns_all' in choice:
        choice = ['burns_muvmass', 'burns_muvz', 'burns_muvmu', 'burns_residualsvz', 'burns_residualsvmass']

    if ('burns_muvmass' in choice) or ('altas_muvmass' in choice):
        try:
            fig, axs = plt.subplots(1, 2, figsize=plot_size, gridspec_kw={'width_ratios': plot_ratio})

            if 'burns_muvmass' in choice:
                fig.suptitle('Host Mass of CSP 91bg-like SNe Ia') # Figure Title
                objs = gen.dict_unpacker(BURNS_SAVE_TXT)
            elif 'altas_muvmass' in choice:
                fig.suptitle('Host Mass of ATLAS 91bg-like SNe Ia') # Figure Title
                objs = gen.dict_unpacker(ATLAS_SAVE_TXT)

            mass_hist = []
            for obj in objs:
                mu, mu_err, mass = float(objs[obj]['mu']), float(objs[obj]['mu_err']), float(objs[obj]['logstellarmass'])
                if mass > 0:
                    axs[0].errorbar(mu, mass, xerr=mu_err, fmt='o')
                    if labels:
                        axs[0].text(mu, mass, obj, size='x-small', va='top')
                    mass_hist.append(mass)
            axs[1].hist(mass_hist, bins=hist_bins[0], orientation="horizontal")

            axs[0].set(xlabel='Distance Modulus', ylabel='Host Mass') # Sub-plot Labels
            axs[1].get_yaxis().set_visible(False)
            plt.tight_layout()
            plt.show()
        except KeyError:
            print("[KeyError] Please run 'ghost_host_galaxy()' before attempting to plot.")

    if ('burns_muvz' in choice) or ('altas_muvz' in choice):
        plot_scale = 'linear'
        sigma = 3

        try:
            fig, axs = plt.subplots(1, 2, figsize=plot_size, gridspec_kw={'width_ratios': plot_ratio})

            if 'burns_muvz' in choice:
                fig.suptitle('Distance Modulus vs. Redshift of CSP 91bg-like SNe Ia)') # Figure Title
                objs = gen.dict_unpacker(BURNS_SAVE_TXT)
            elif 'altas_muvz' in choice:
                fig.suptitle('Distance Modulus vs. Redshift of ATLAS 91bg-like SNe Ia') # Figure Title
                objs = gen.dict_unpacker(ATLAS_SAVE_TXT)

            # Plot -- mu vs z
            mu_hist = []
            for obj in objs:
                mu, mu_err, z = float(objs[obj]['mu']), float(objs[obj]['mu_err']), float(objs[obj]['z'])
                if not raw and mu <= 0:
                    continue
                axs[0].errorbar(z, mu, yerr=mu_err*sigma, fmt='o')
                if labels:
                    axs[0].text(z, mu, obj, size='x-small', va='top')
                mu_hist.append(mu)
            axs[1].hist(mu_hist, bins=hist_bins[1], orientation="horizontal")

            # Plot cosmology -- Use the distmod method of the cosmology object
            z_cosmo = np.linspace(0.001, 0.08, 100)
            axs[0].plot(z_cosmo, COSMO_MODEL.distmod(z_cosmo), alpha=0.4, linewidth=5,
                        label='H0:'+str(COSMO_MODEL._H0)+'\nOm0: '+str(COSMO_MODEL._Om0)+'\nOde0: '+str(COSMO_MODEL._Ode0))

            # Formating
            axs[0].set(xlabel='Redshift', ylabel='Distance Modulus') # Sub-plot Labels
            axs[0].set_xscale(plot_scale); axs[0].set_yscale(plot_scale)
            axs[0].legend()
            axs[1].get_yaxis().set_visible(False)
            fig.tight_layout()

            plt.show()
        except KeyError:
            print("[KeyError] Please run 'snpy_fit()' before attempting to plot.")

    if ('burns_muvmu' in choice) or ('altas_muvmu' in choice):
        plot_scale = 'linear'
        plt.figure(figsize=plot_size)

        try:
            if 'burns_muvmu' in choice:
                plt.title('SNooPy Distance Modulus vs. Cosmological Distance Modulus of CSP 91bg-like SNe Ia)') # Figure Title
                objs = gen.dict_unpacker(BURNS_SAVE_TXT)
            elif 'altas_muvmu' in choice:
                plt.title('SNooPy Distance Modulus vs. Cosmological Distance Modulus of ATLAS 91bg-like SNe Ia') # Figure Title
                objs = gen.dict_unpacker(ATLAS_SAVE_TXT)

            # Plot -- mu vs mu
            x_mu, y_mu = [], []
            for obj in objs:
                mu_snpy, z = float(objs[obj]['mu']), float(objs[obj]['z'])
                mu_cosmo = COSMO_MODEL.distmod(z).value
                plt.scatter(mu_snpy, mu_cosmo, marker='o')
                x_mu.append(mu_snpy)
                y_mu.append(mu_cosmo)
                if labels:
                    axs[0].text(mu_snpy, mu_cosmo, obj, size='x-small', va='top')

            # Plot trendline
            z = np.polyfit(x_mu, y_mu, 1)
            p = np.poly1d(z)
            plt.plot(x_mu, p(x_mu), alpha=0.4, linewidth=5, label='Trendline: '+"y=%.6fx+(%.6f)"%(z[0],z[1]))

            # Formating
            plt.xlabel('SNooPy Distance Modulus'); plt.ylabel('Cosmological Distance Modulus')
            plt.legend()

            plt.show()
        except KeyError:
            print("[KeyError] Please run 'snpy_fit()' before attempting to plot.")

    if ('burns_residualsvz' in choice) or ('altas_residualsvz' in choice):
        sigma = 0

        try:
            fig, axs = plt.subplots(1, 2, figsize=plot_size, gridspec_kw={'width_ratios': plot_ratio})
            if 'burns_residualsvz' in choice:
                fig.suptitle('Hubble Residuals vs. Redshift of CSP 91bg-like SNe Ia\n '+
                             'Dist. Sigma: '+str(sigma)) # Figure Title
                objs = gen.dict_unpacker(BURNS_SAVE_TXT)
            elif 'altas_residualsvz' in choice:
                fig.suptitle('Hubble Residuals vs. Redshift of ATLAS 91bg-like SNe Ia\n'+
                             'Dist. Sigma: '+str(sigma)) # Figure Title
                objs = gen.dict_unpacker(ATLAS_SAVE_TXT)

            # Compute mu_cosmo
            mu_cosmo, mu_snpy, mu_err, z, objnames = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
            for obj in objs:
                objnames = np.append(objnames, obj)
                z = np.append(z, float(objs[obj]['z']))
                mu_err = np.append(mu_err, float(objs[obj]['mu_err']))
                mu_snpy = np.append(mu_snpy, float(objs[obj]['mu']))
                mu_cosmo = np.append(mu_cosmo, COSMO_MODEL.distmod(float(objs[obj]['z'])).value)
            mu_res = (mu_snpy - mu_cosmo) - np.average(mu_snpy - mu_cosmo)

            # Plot -- mu vs z
            mu_hist = []
            for n in range(len(objnames)):
                if not raw and z[n] < 0.01:
                    continue
                axs[0].errorbar(z[n], mu_res[n], yerr=mu_err[n]*sigma, label=objnames[n], fmt="o")
                mu_hist.append(mu_res[n])
                if labels:
                    axs[0].text(z[n], mu_res[n], objnames[n], size='x-small', va='top')
            axs[1].hist(mu_hist, bins=hist_bins[2], orientation="horizontal")

            # Formatting
            axs[0].set(xlabel='Redshift', ylabel='Hubble Residuals') # Sub-plot Labels
            axs[1].get_yaxis().set_visible(False)
            plt.tight_layout()

            # Limits
            ylimiter = (np.max(np.abs(mu_hist)) + np.max(mu_err*sigma)) + 0.01
            axs[0].set_ylim(-ylimiter, ylimiter); axs[1].set_ylim(-ylimiter, ylimiter)

            plt.show()
        except KeyError:
            print("[KeyError] Please run 'snpy_fit()' before attempting to plot.")

    if ('burns_residualsvmass' in choice) or ('altas_residualsvmass' in choice):
        sigmas = [1, 1]
        objsRemove = ['2022skw', '2023cvq', '2021cad', '2020nta']

        try:
            fig, axs = plt.subplots(1, 2, figsize=plot_size, gridspec_kw={'width_ratios': plot_ratio})

            if 'burns_residualsvmass' in choice:
                fig.suptitle('Hubble Residuals vs. Host Mass of CSP 91bg-like SNe Ia\n '+
                             'Dist. Sigma: '+str(sigmas[0])+' | Mass Sigma: '+str(sigmas[1])) # Figure Title
                objs = gen.dict_unpacker(BURNS_SAVE_TXT)
            elif 'altas_residualsvmass' in choice:
                fig.suptitle('Hubble Residuals vs. Host Mass of ATLAS 91bg-like SNe Ia\n '+
                             'Dist. Sigma: '+str(sigmas[0])+' | Mass Sigma: '+str(sigmas[1])) # Figure Title
                objs = gen.dict_unpacker(ATLAS_SAVE_TXT)

            # Compute mu_cosmo
            mu_cosmo, mu_snpy, mu_err, mass, mass_err, objnames = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
            for obj in objs:
                if obj in objsRemove:
                    continue
                objnames = np.append(objnames, obj)
                mass = np.append(mass, float(objs[obj]['logstellarmass']))
                mass_err = np.append(mass_err, float(objs[obj]['logstellarmasserr']))
                mu_err = np.append(mu_err, float(objs[obj]['mu_err']))
                mu_snpy = np.append(mu_snpy, float(objs[obj]['mu']))
                mu_cosmo = np.append(mu_cosmo, COSMO_MODEL.distmod(float(objs[obj]['z'])).value)
            mu_res = (mu_snpy - mu_cosmo) - np.average(mu_snpy - mu_cosmo)

            # Plot -- mu vs z
            mu_hist = []
            for n in range(len(objnames)):
                if not raw and mass[n] <= 0:
                    continue
                axs[0].errorbar(mass[n], mu_res[n], xerr=mass_err[n]*sigmas[1], yerr=mu_err[n]*sigmas[0], fmt='o')
                mu_hist.append(mu_res[n])
                if labels:
                    axs[0].text(mass[n], mu_res[n], objnames[n], size='x-small', va='top')
            axs[1].hist(mu_hist, bins=hist_bins[3], orientation="horizontal")

            # Formatting
            axs[0].set(xlabel='Host Mass', ylabel='Hubble Residuals') # Sub-plot Labels
            axs[1].get_yaxis().set_visible(False)
            # axs[0].invert_xaxis()
            plt.tight_layout()

            # Limits
            ylimiter = (np.max(np.abs(mu_hist)) + np.max(mu_err*sigmas[0]))
            ylimiter = 1.5
            axs[0].set_ylim(-ylimiter, ylimiter); axs[1].set_ylim(-ylimiter, ylimiter)
            # axs[0].set_xlim(9.75, 11.5); axs[1].set_xlim(9.75, 11.5)

            plt.show()
        except KeyError:
            print("[KeyError] Please run 'snpy_fit()' before attempting to plot.")

    if ('burns_res_zcorr' in choice) or ('altas_res_zcorr' in choice):
        sigma = 0
        objsRemove = ['2022skw', '2023cvq']

        try:
            fig, axs = plt.subplots(1, 2, figsize=plot_size, gridspec_kw={'width_ratios': plot_ratio})
            if 'burns_res_zcorr' in choice:
                fig.suptitle('Hubble Residuals vs. Redshift of CSP 91bg-like SNe Ia\n '+
                             'Dist. Sigma: '+str(sigma)) # Figure Title
                objs = gen.dict_unpacker(BURNS_SAVE_TXT)
            elif 'altas_res_zcorr' in choice:
                fig.suptitle('Hubble Residuals vs. Redshift of ATLAS 91bg-like SNe Ia\n'+
                             'Dist. Sigma: '+str(sigma)) # Figure Title
                objs = gen.dict_unpacker(ATLAS_SAVE_TXT)

            # Compute mu_cosmo
            mu_cosmo, mu_snpy, mu_err, z, objnames = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
            for obj in objs:
                if not raw and float(objs[obj]['z']) < 0.01:
                    continue
                elif obj in objsRemove:
                    continue
                objnames = np.append(objnames, obj)
                z = np.append(z, float(objs[obj]['z']))
                mu_err = np.append(mu_err, float(objs[obj]['mu_err']))
                mu_snpy = np.append(mu_snpy, float(objs[obj]['mu']))
                mu_cosmo = np.append(mu_cosmo, COSMO_MODEL.distmod(float(objs[obj]['z'])).value)
            mu_res = (mu_snpy - mu_cosmo) - np.average(mu_snpy - mu_cosmo)

            # Plot -- mu vs z
            mu_hist = []
            for n in range(len(objnames)):

                axs[0].errorbar(z[n], mu_res[n], yerr=mu_err[n]*sigma, label=objnames[n], fmt="o")
                mu_hist.append(mu_res[n])
                if labels:
                    axs[0].text(z[n], mu_res[n], objnames[n], size='x-small', va='top')
            axs[1].hist(mu_hist, bins=hist_bins[2], orientation="horizontal")

            # Formatting
            axs[0].set(xlabel='Redshift', ylabel='Hubble Residuals') # Sub-plot Labels
            axs[1].get_yaxis().set_visible(False)
            plt.tight_layout()

            # Limits
            ylimiter = (np.max(np.abs(mu_hist)) + np.max(mu_err*sigma)) + 0.1
            axs[0].set_ylim(-ylimiter, ylimiter); axs[1].set_ylim(-ylimiter, ylimiter)

            # Redshift correction
            # chi2 = malmquist_bias_corr(mu_res, z, mu_err, m, b)






            plt.close()
            # plt.show()
        except KeyError:
            print("[KeyError] Please run 'snpy_fit()' before attempting to plot.")


    return
