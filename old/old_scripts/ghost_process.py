from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from astropy import cosmology as cosmo
from astro_ghost.ghostHelperFunctions import getTransientHosts
import shutil
import numpy as np

from astroquery.sdss import SDSS
from astropy import coordinates as coords
from astroquery.mast import Catalogs


import matplotlib.pyplot as plt

import glob

import scripts.general as gen
import scripts.atlas_process as atlas
import scripts.ztf_process as ztf

COSMO_MODEL = cosmo.FlatLambdaCDM(H0=70, Om0=0.3)

GHOST_DATA = '../data/GHOST/'
TEST_ROOT = '../tests/'

ATLAS_SAVE_TXT = '../snpy/atlas/atlas_saved.txt'
BURNS_SAVE_TXT = '../snpy/burns/burns_saved.txt'
ZTF_SAVE_TXT = '../snpy/ztf/ztf_saved.txt'
COMBINED_SAVE_TXT = '../snpy/combined_saved.txt'
DR3_SAVE_TXT = '../txts/DR3_fits.dat'

TNS_KEY_TXT = '../working_data/TNS_key.txt'

def ghost_host_galaxy(dict_path, save_loc=TEST_ROOT, keep_data=True, update_saved=False):
    cosmo = FlatLambdaCDM(70, 0.3) # Hubble Constant, Omega-Matter
    data = gen.dict_unpacker(dict_path)
    all_z, all_logstellarmass = [], []
    print('[+++] Finding host galaxy mass using GHOST...')

    i = 0
    for obj in data:

        # i += 1
        # if i < 1:
        #     continue
        # elif i > 1:
        #     break


        ra, dec, z = float(data[obj]['ra']), float(data[obj]['dec']), float(data[obj]['z'])

        print('\n[', list(data).index(obj)+1, '/', len(data), ']', obj, '|', data[obj]['ra'], data[obj]['dec'], data[obj]['z'])
        print('---------------------------------------------------------------------------------------------------------')

        transient_position = SkyCoord(ra, dec, unit=u.deg)
        try:
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                host_data = getTransientHosts(transientCoord=[transient_position], transientName=[obj], verbose=False, starcut="gentle", savepath=save_loc+'ghost_stuff/', GHOSTpath=GHOST_DATA)

            # If GLADE can't get magnitudes -- using NED name to get SDSS data
            print('Identified Host Galaxy:', host_data.loc[0, 'NED_name'])
            if np.isnan(host_data.loc[0, 'gKronMag']) or np.isnan(host_data.loc[0, 'iKronMag']):
                print('GLADE does not contain the g-mag & i-mag, searching SDSS...')
                host_position = SkyCoord(host_data['raMean'], host_data['decMean'], unit=u.deg)
                result = SDSS.query_crossid(host_position, photoobj_fields=['modelMag_g', 'modelMagErr_g', 'modelMag_i', 'modelMagErr_i'])

                # If GLADE or SDSS has magnitudes -- using PanSTARRS
                if result == None:
                    print('GLADE & SDSS do not contain the g-mag & i-mag, searching PanSTARRS...')
                    catalog_data = Catalogs.query_region(transient_position, radius=2 * u.arcsec, catalog="Panstarrs")
                    # print(catalog_data['rMeanKronMag'], catalog_data['iMeanKronMag'])

                    if len(catalog_data) == 0 or isinstance(catalog_data['gMeanKronMagErr'].value[0], np.ma.core.MaskedConstant) or isinstance(catalog_data['iMeanKronMagErr'].value[0], np.ma.core.MaskedConstant):
                        print('[!!!!!] GLADE, SDSS, and PanSTARRS do not contain the g-mag & i-mag data...')
                        gMag, iMag, iAbsMag = np.nan, np.nan, np.nan
                        gMagErr, iMagErr, iAbsMagErr = np.nan, np.nan, np.nan
                    else:
                        gMag, iMag, iAbsMag = catalog_data['gMeanKronMag'].value[0], catalog_data['gMeanKronMag'].value[0], (catalog_data['iMeanKronMag'].value[0] - cosmo.distmod(z).value)
                        gMagErr, iMagErr, iAbsMagErr = catalog_data['gMeanKronMagErr'].value[0], catalog_data['iMeanKronMagErr'].value[0], catalog_data['iMeanKronMagErr'].value[0]
                else:
                    # SDSS Results
                    gMag, iMag, iAbsMag = result['modelMag_g'].value[0], result['modelMag_i'].value[0], (result['modelMag_i'].value[0] - cosmo.distmod(z).value)
                    gMagErr, iMagErr, iAbsMagErr = result['modelMagErr_g'].value[0], result['modelMagErr_i'].value[0], result['modelMagErr_i'].value[0]
            else:
                # GLADE results
                gMag, iMag, iAbsMag = host_data['gKronMag'].loc[0], host_data['iKronMag'].loc[0], (host_data['iKronMag'].loc[0] - cosmo.distmod(z).value)
                gMagErr, iMagErr, iAbsMagErr = host_data['gKronMagErr'].loc[0], host_data['iKronMagErr'].loc[0], host_data['iKronMagErr'].loc[0]




            #  Mass Calculation -- Taylor et. al. 2011 -- eq. 8
            logstellarmass = (1.15 + (0.7*(gMag - iMag)) - (0.4*(iAbsMag)))

            # Error Propogation
            giMagErr = np.sqrt((gMagErr ** 2) + (iMagErr ** 2))
            logstellarmasserr = np.sqrt(((0.7**2)*(giMagErr**2)) + ((0.4**2)*(iAbsMagErr**2)))


        except Exception as error:
            print(obj+':', error)
            logstellarmass, logstellarmasserr = 0.00, 0.00

        if np.isnan(logstellarmass) == False and logstellarmass > 0.00:
            print('Success!', obj, 'host galaxy has a mass of:', logstellarmass, '+/-', logstellarmasserr, 'logM_* / [Msun]')
            all_z.append(z)
            all_logstellarmass.append(logstellarmass)
            if update_saved:
                data[obj].update({'logstellarmass': logstellarmass, 'logstellarmasserr': logstellarmasserr})
        else:
            print('Failed to find host galaxy!')
            if update_saved:
                data[obj].update({'logstellarmass': 0.00, 'logstellarmasserr': 0.00})

        # if i > 0:
        #     break
        # else:
        #     i =+ 1


    print('\nSuccessfully found mass of', len(all_z), '/', len(data), 'host galaxies!')
    if not keep_data:
        print('Removing GHOST data...')
        shutil.rmtree(save_loc+'ghost_stuff/') # Clear messy data
    if update_saved:
        print('Saving data to'+dict_path+'...')
        gen.dict_packer(data, dict_path)

    return all_z, all_logstellarmass
def combined_data():
    atlas_objs = gen.dict_unpacker(ATLAS_SAVE_TXT)
    burns_objs = gen.dict_unpacker(BURNS_SAVE_TXT)
    ztf_objs = gen.dict_unpacker(ZTF_SAVE_TXT)

    combined_objs = {}
    for obj in burns_objs:
        combined_objs.update({obj: burns_objs[obj]})
    for obj in atlas_objs:
        if obj not in combined_objs:
            combined_objs.update({obj: atlas_objs[obj]})

    gen.dict_packer(combined_objs, '../snpy/combined_saved.txt')

    return

def atlas_ztf_joint_fitting():
    ATLASobjs = atlas.atlas_processing(err_max=1000, n_iter=0, sleep_t=5, use_TNS=True, loc_TNS=TNS_KEY_TXT)
    ZTFobjs = ztf.ztf_alt_processing(paths=glob.glob('../data/ZTF/*.txt'), min_pts=10, mag_err_max=0.75, flux_err_max=80)

    print(ZTFobjs['2020nta'].keys())

    # combined_objs = ATLASobjs
    # for obj in ZTFobjs:
    #     if obj in combined_objs:
    #         print(obj)


    return

def ghost_plotting(choice, plot_size = (18, 6), plot_ratio = [10, 1], hist_bins = [50, 50, 10],
                   labels = True, raw = False, save = False, ignore_type=None):
    badObjs = ['2018ame', '2018ast', '2020nta', '2021agej', '2021agnf', '2021cad', '2021mab', '2022an', '2022fjx',
               '2022omn', '2022oux', '2022skw', '2023cvq', '2023ex', '2023mkp', '2023omo', '2023yrs']
    okayObjs = ['2018jag', '2019cp', '2021bls', '2021gel', '2021jbp', '2019be', '2020acoo', '2021jvp', '2021oyx',
                '2021pom',
                '2021uve', '2022aaok', '2022abom', '2022dsu', '2022rjs', '2022ydr', '2022yv', '2023abdv', '2023fwb',
                '2023jah']
    goodObjs = ['2019exc', '2019cdc', '2019ecx', '2020fhs', '2021fnr', '2019ecx', '2022aecb', '2022ihz', '2022xhh',
                '2022xkq', '2023acdv', '2023bhm', '2024bjb', '2024jhk']

    if 'atlas_all' in choice:
        choice = ['altas_muvmass', 'altas_muvz', 'altas_muvmu', 'altas_residualsvz', 'altas_residualsvmass']
    if 'burns_all' in choice:
        choice = ['burns_muvmass', 'burns_muvz', 'burns_muvmu', 'burns_residualsvz', 'burns_residualsvmass']

    if 'burns_dr3' in choice:
        burns_objs = gen.dict_unpacker(BURNS_SAVE_TXT)
        dr3_objs = gen.dict_unpacker(DR3_SAVE_TXT, delimiter=None)

        fig, axs = plt.subplots(1, 3, figsize=(18, 6), gridspec_kw={'width_ratios': [1, 1, 1]})

        for obj in burns_objs:
            axs[0].errorbar(float(burns_objs[obj]['st']), float(dr3_objs['SN'+obj]['st']),
                            xerr=float(burns_objs[obj]['st_err']), yerr=float(dr3_objs['SN'+obj]['e_st']),
                            fmt='o', linewidth=0.5)
            axs[1].errorbar(float(burns_objs[obj]['Tmax'])+53000, float(dr3_objs['SN'+obj]['Tmax']),
                            xerr=float(burns_objs[obj]['Tmax_err']), yerr=float(dr3_objs['SN'+obj]['eTmax']),
                            fmt='o', linewidth=0.5)
            axs[2].errorbar(float(burns_objs[obj]['EBVhost']), float(dr3_objs['SN'+obj]['EBVhost']),
                            xerr=float(burns_objs[obj]['EBVhost_err']), yerr=float(dr3_objs['SN'+obj]['e_EBVhost']),
                            fmt='o', linewidth=0.5)

            # Labels
            if labels:
                axs[0].text(float(burns_objs[obj]['st']), float(dr3_objs['SN'+obj]['st']), obj, size='x-small', va='top')
                axs[1].text(float(burns_objs[obj]['Tmax'])+53000, float(dr3_objs['SN'+obj]['Tmax']), obj, size='x-small', va='top')
                axs[2].text(float(burns_objs[obj]['EBVhost']), float(dr3_objs['SN'+obj]['EBVhost']), obj, size='x-small', va='top')

            # # Output
            # print('ST: ', float(burns_objs[obj]['st']), '+/-', float(burns_objs[obj]['st_err']), '|',
            #       float(dr3_objs['SN'+obj]['st']), '+/-', float(dr3_objs['SN'+obj]['e_st']))
            # print('Tmax: ', float(burns_objs[obj]['Tmax'])+53000, '+/-', float(burns_objs[obj]['Tmax_err']), '|',
            #       float(dr3_objs['SN'+obj]['Tmax']), '+/-', float(dr3_objs['SN'+obj]['eTmax']))
            # print('EBVhost: ', float(burns_objs[obj]['EBVhost']), '+/-', float(burns_objs[obj]['EBVhost_err']), '|',
            #       float(dr3_objs['SN'+obj]['EBVhost']), '+/-', float(dr3_objs['SN'+obj]['e_EBVhost']))

        # 1-to-1 lines
        axs[0].axline((0.2,0.2), slope=1, alpha=0.3, linewidth=4)
        axs[1].axline((53400,53400), slope=1, alpha=0.3, linewidth=4)
        axs[2].axline((0,0), slope=1, alpha=0.3, linewidth=4)

        axs[0].set_title('Stretch [st]'); axs[1].set_title('SNooPy vs. Burns Parameters\n Tmax [MJD+53000]'); axs[2].set_title('EBVhost')
        axs[0].set(ylabel='Burns'); axs[1].set(xlabel='SNooPy')
        plt.tight_layout()

        if save:
            plt.savefig('../save/plots/SNooPy_v_Burns.png')

        plt.show()

    if ('burns_muvmass' in choice) or ('atlas_muvmass' in choice):
        try:
            fig, axs = plt.subplots(1, 2, figsize=plot_size, gridspec_kw={'width_ratios': plot_ratio})

            if 'burns_muvmass' in choice:
                title = 'Host Mass of CSP 91bg-like SNe Ia'
                objs = gen.dict_unpacker(BURNS_SAVE_TXT)
            elif 'atlas_muvmass' in choice:
                title = 'Host Mass of ATLAS 91bg-like SNe Ia'
                objs = gen.dict_unpacker(ATLAS_SAVE_TXT)
            fig.suptitle(title)  # Figure Title

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

            if save:
                fig.savefig('../save/plots/dist_v_mass-'+choice[0][:5]+'.png')

            plt.show()

        except KeyError:
            print("[KeyError] Please run 'ghost_host_galaxy()' before attempting to plot.")

    if ('burns_muvz' in choice) or ('atlas_muvz' in choice):
        plot_scale = 'linear'
        sigma = 3

        try:
            fig, axs = plt.subplots(1, 2, figsize=plot_size, gridspec_kw={'width_ratios': plot_ratio})

            if 'burns_muvz' in choice:
                fig.suptitle('Distance Modulus vs. Redshift of CSP 91bg-like SNe Ia)') # Figure Title
                objs = gen.dict_unpacker(BURNS_SAVE_TXT)
            elif 'atlas_muvz' in choice:
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

            if save:
                fig.savefig('../save/plots/dist_v_z-'+choice[0][:5]+'.png')

            plt.show()

        except KeyError:
            print("[KeyError] Please run 'snpy_fit()' before attempting to plot.")

    if ('burns_muvmu' in choice) or ('atlas_muvmu' in choice):
        plot_scale = 'linear'
        plt.figure(figsize=plot_size)

        try:
            if 'burns_muvmu' in choice:
                plt.title('SNooPy Distance Modulus vs. Cosmological Distance Modulus of CSP 91bg-like SNe Ia)') # Figure Title
                objs = gen.dict_unpacker(BURNS_SAVE_TXT)
            elif 'atlas_muvmu' in choice:
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

            if save:
                fig.savefig('../save/plots/dist_v_dist-'+choice[0][:5]+'.png')

            plt.show()

        except KeyError:
            print("[KeyError] Please run 'snpy_fit()' before attempting to plot.")

    if ('burns_residualsvz' in choice) or ('atlas_residualsvz' in choice) or ('combined_residualsvz' in choice):
        sigma = 1
        objsRemove = ['2021cad', '2020nta', '2022skw', '2023cvq', '2023yrs', '2018ast']

        try:
            fig, axs = plt.subplots(1, 2, figsize=plot_size, gridspec_kw={'width_ratios': plot_ratio})
            if 'burns_residualsvz' in choice:
                fig.suptitle('Hubble Residuals vs. Redshift of CSP 91bg-like SNe Ia\n '+
                             'Dist. Sigma: '+str(sigma)) # Figure Title
                objs = gen.dict_unpacker(BURNS_SAVE_TXT)
            elif 'atlas_residualsvz' in choice:
                fig.suptitle('Hubble Residuals vs. Redshift of ATLAS 91bg-like SNe Ia\n'+
                             'Dist. Sigma: '+str(sigma)) # Figure Title
                objs = gen.dict_unpacker(ATLAS_SAVE_TXT)
            elif 'combined_residualsvz' in choice:
                fig.suptitle('Hubble Residuals vs. Redshift of ATLAS+CSP 91bg-like SNe Ia\n'+
                             'Dist. Sigma: '+str(sigma)) # Figure Title
                objs = gen.dict_unpacker(COMBINED_SAVE_TXT)

            # Compute mu_cosmo
            mu_cosmo, mu_snpy, mu_err, z, objnames = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
            for obj in objs:
                if ('bad' in ignore_type) and (obj in badObjs):
                    continue
                elif ('okay' in ignore_type) and (obj in okayObjs):
                    continue
                elif ('good' in ignore_type) and (obj in goodObjs):
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

            # Scatter
            print('Scatter: ', np.std(mu_res))

            if save:
                fig.savefig('../save/plots/resid_v_z-'+choice[0][:5]+'.png')

            plt.show()

        except KeyError:
            print("[KeyError] Please run 'snpy_fit()' before attempting to plot.")

    if ('burns_residualsvmass' in choice) or ('atlas_residualsvmass' in choice) or ('combined_residualsvmass' in choice):
        sigmas = [1, 1]

        try:
            fig, axs = plt.subplots(1, 2, figsize=plot_size, gridspec_kw={'width_ratios': plot_ratio})

            if 'burns_residualsvmass' in choice:
                fig.suptitle('Hubble Residuals vs. Host Mass of CSP 91bg-like SNe Ia\n '+
                             'Dist. Sigma: '+str(sigmas[0])+' | Mass Sigma: '+str(sigmas[1])) # Figure Title
                objs = gen.dict_unpacker(BURNS_SAVE_TXT)
            elif 'atlas_residualsvmass' in choice:
                fig.suptitle('Hubble Residuals vs. Host Mass of ATLAS 91bg-like SNe Ia\n '+
                             'Dist. Sigma: '+str(sigmas[0])+' | Mass Sigma: '+str(sigmas[1])) # Figure Title
                objs = gen.dict_unpacker(ATLAS_SAVE_TXT)
            elif 'combined_residualsvmass' in choice:
                fig.suptitle('Hubble Residuals vs. Host Mass of ATLAS+CSP 91bg-like SNe Ia\n '+
                             'Dist. Sigma: '+str(sigmas[0])+' | Mass Sigma: '+str(sigmas[1])) # Figure Title
                objs = gen.dict_unpacker(COMBINED_SAVE_TXT)

            # Compute mu_cosmo
            mu_cosmo, mu_snpy, mu_err, mass, mass_err, objnames = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
            for obj in objs:
                if ('bad' in ignore_type) and (obj in badObjs):
                    continue
                elif ('okay' in ignore_type) and (obj in okayObjs):
                    continue
                elif ('good' in ignore_type) and (obj in goodObjs):
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
            axs[0].set_ylim(-ylimiter, ylimiter); axs[1].set_ylim(-ylimiter, ylimiter)

            # Scatter
            print('Scatter: ', np.std(mu_res))

            if save:
                fig.savefig('../save/plots/resid_v_mass-' + choice[0][:5] + '.png')

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
