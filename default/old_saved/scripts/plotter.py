import os
import sys
import snpy
import time as systime
import numpy as np
import matplotlib.pyplot as plt

import scripts.general as gen

def lc_plot(objs, y_type = 'flux', pause_time=2, color_wheel = ['orange', 'cyan', 'violet', 'red', 'blue'],
            quiet=False, save_plots=True, save_loc='../snpy/misc_plots/'):
    print('[+++] Plotting LC data...')
    color_wheel = [None, None, None, None, None, None, None, None, None, None, None]

    # Check quiet
    if quiet:
        sys.stdout = open(os.devnull, 'w')

    for obj in objs:
        print('-------------------------------------------------------------------------------------------------------')
        print('[', list(objs).index(obj)+1, '/', len(objs), '] -', obj)

        plt.figure(figsize=(12, 6))

        # Get list of all filters
        filter_wheel = []
        for f_w in objs[obj]['filters']:
            if f_w not in filter_wheel:
                filter_wheel.append(f_w)

        for f_w in filter_wheel:
            f_indexs = np.where(objs[obj]['filters'] == f_w)[0]
            plt.errorbar(objs[obj]['time'][f_indexs], objs[obj][y_type][f_indexs], yerr=objs[obj]['d'+y_type][f_indexs],
                        color=color_wheel[filter_wheel.index(f_w)], label=f_w, fmt='o')
        systime.sleep(pause_time)

        if y_type == 'mag':
            plt.gca().invert_yaxis()
        plt.title(obj)
        plt.xlabel('JD'); plt.ylabel(y_type)
        plt.legend()
        if save_plots:
            plt.savefig(save_loc + obj + '_lc.png')
            print(obj, '-- Plot saved to', save_loc + obj + '_lc.png')
        plt.show()
        plt.close()


    # Restore print statements
    sys.stdout = sys.__stdout__

    return
def lc_replot(lc_path, save_plot=False, save_loc='../default/', colors=None, spread=None, stacked=True):
    n_s = snpy.get_sn(lc_path)
    # if colors is None:
    #     plt_args = {'single': stacked, 'colors': colors}
    # else:
    plt_args = {'single': stacked}

    if spread != None:
        plt_args.update({'xrange': (n_s.parameters['Tmax']-spread[0], n_s.parameters['Tmax']+spread[1])})
    if save_plot:
        plt_args.update({'outfile': save_loc + n_s.name + '_snpylc.png'})
    n_s.plot(**plt_args)
    plt.show()
    return
def residual_plotter(path, x_params, sigma=[1,1], labels=False, raw=False, extra_info=False, ignore_type=[], save_plot=False, save_loc='../default/'):
    # Get reviewed fits
    with open(gen.get_constants()['reviewed_fits_txt'], 'r') as f:
        reviewed_good_fits, reviewed_okay_fits, reviewed_bad_fits = [], [], []
        for container in [reviewed_good_fits, reviewed_okay_fits, reviewed_bad_fits]:
            for i in f.readline().split(', ')[1:]:
                container.append(i[1:-1])
            container[-1] = container[-1][:-1]

    # Pull data from saved text
    objs = gen.dict_handler(choice='unpack', path=path)

    fig, axs = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'width_ratios': [10, 1]}, constrained_layout=True)
    atlas_key, ztf_key, csp_key = True, True, True
    colors = {'ATLAS': 'blue', 'ZTF': 'red', 'CSP': 'orange'}
    mu_res_hist = []
    for obj in objs:
        if x_params[0] == 'host_mass':
            n_x_err = float(objs[obj][x_params[0]+'_err'])
        else:
            n_x_err = 0.00
        n_x = float(objs[obj][x_params[0]])
        n_mu_res = float(objs[obj]['mu']) - gen.current_cosmo().distmod(float(objs[obj]['z'])).value
        n_mu_err = float(objs[obj]['mu_err'])

        # Clean data
        if not raw:
            objname = obj
            if objname[:2] == 'SN':
                objname = objname[2:]
            if ('good' in ignore_type) and (objname in reviewed_good_fits):
                continue
            if ('okay' in ignore_type) and (objname in reviewed_okay_fits):
                continue
            if ('bad' in ignore_type) and (objname in reviewed_bad_fits):
                continue
            if x_params[0] == 'host_mass' and n_x == 0:  # Remove null masses
                continue

        # Save for histogram
        mu_res_hist.append(n_mu_res)

        # Plot points
        if path.split('/')[-1][:-10] == 'combined':
            if atlas_key and (objs[obj]['origin'] == 'ATLAS'):
                axs[0].errorbar(n_x, n_mu_res, yerr=n_mu_err * sigma[0], xerr=n_x_err * sigma[1],
                                fmt='o', color=colors[objs[obj]['origin']], label=objs[obj]['origin'])
                atlas_key = False
            elif ztf_key and (objs[obj]['origin'] == 'ZTF'):
                axs[0].errorbar(n_x, n_mu_res, yerr=n_mu_err * sigma[0], xerr=n_x_err * sigma[1],
                                fmt='o', color=colors[objs[obj]['origin']], label=objs[obj]['origin'])
                ztf_key = False
            elif csp_key and (objs[obj]['origin'] == 'CSP'):
                axs[0].errorbar(n_x, n_mu_res, yerr=n_mu_err * sigma[0], xerr=n_x_err * sigma[1],
                                fmt='o', color=colors[objs[obj]['origin']], label=objs[obj]['origin'])
                csp_key = False
            else:
                axs[0].errorbar(n_x, n_mu_res, yerr=n_mu_err * sigma[0], xerr=n_x_err * sigma[1],
                                fmt='o', color=colors[objs[obj]['origin']])
        else:
            axs[0].errorbar(n_x, n_mu_res, yerr=n_mu_err * sigma[0], xerr=n_x_err * sigma[1],
                            fmt='o')

        # Labels
        if labels:
            axs[0].text(n_x, n_mu_res, obj, size='x-small', va='top')

    # Histogram
    axs[1].hist(mu_res_hist, bins=40, orientation="horizontal")

    # Formatting
    ylimiter = np.max(np.abs(mu_res_hist))+0.5
    axs[0].set_ylim(-ylimiter, ylimiter); axs[1].set_ylim(-ylimiter, ylimiter)

    if extra_info:
        fig.suptitle("Hubble Residuals vs. " + x_params[1] + " of '"+(path.split('/')[-1].split('_')[0]).upper()+"' 91bg-like SNe Ia\n" +  # Figure Title
                     'Dist. Sigma: ' + str(sigma[0]) + ' | ' + x_params[1] + ' Sigma: ' + str(sigma[0]) +
                     ' | Scatter: ' + str(round(np.std(mu_res_hist), 2)) + ' | # of pts: ' + str(len(mu_res_hist)), size='medium')
    else:
        fig.suptitle("Hubble Residuals vs. " + x_params[1] + " of '"+(path.split('/')[-1][:-10]).upper()+"' 91bg-like SNe Ia")
    axs[0].set(xlabel=x_params[1], ylabel='Hubble Residuals')  # Sub-plot Labels
    axs[1].get_yaxis().set_visible(False) # Turn off y-axis labels
    if path.split('/')[-1][:-10] == 'combined':
        axs[0].legend()
    if save_plot:
        plt.savefig(save_loc+'hubble_res_v_'+x_params[1]+'.png')
    plt.show()
    return
def alt_residual_plotter(path, x_params=['z', 'Redshift']):
    # Pull data from saved text
    objs = gen.dict_handler(choice='unpack', path=path)

    data_set = path.split('/')[-1].split('_')[0].upper()

    fig, axs = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'width_ratios': [10, 1]}, constrained_layout=True)
    mu_hist_arr = np.array([])
    color_wheel = {'ZTF': 'red', 'ATLAS': 'blue', 'CSP': 'orange'}
    key = {'ZTF': True, 'ATLAS': True, 'CSP': True}
    for obj in objs:
        objname = obj

        n_x, n_x_err = float(objs[obj][x_params[0]]), 0.00
        n_y = float(objs[obj]['mu']) - gen.current_cosmo().distmod(float(objs[obj]['z'])).value
        n_y_err = float(objs[obj]['mu_err'])
        color = None
        if x_params[0] == 'host_mass':
            n_x_err = float(objs[obj][x_params[0] + '_err'])
        if data_set == 'COMBINED':
            color = color_wheel[objs[obj]['origin']]


        # Plot points
        axs[0].errorbar(n_x, n_y, xerr=n_x_err, yerr=n_y_err, color=color, fmt='o')
        mu_hist_arr = np.append(mu_hist_arr, n_y)

    # Plot histogram
    axs[1].hist(mu_hist_arr, bins=40, orientation="horizontal")

    # Formatting
    ylimiter = np.max(np.abs(mu_hist_arr))+0.5
    axs[0].set_ylim(-ylimiter, ylimiter); axs[1].set_ylim(-ylimiter, ylimiter)
    fig.suptitle("Hubble Residuals vs. " + x_params[1] + " of 'ATLAS' 91bg-like SNe Ia")
    axs[0].set(xlabel="Redshift", ylabel='Hubble Residuals')  # Sub-plot Labels
    axs[1].get_yaxis().set_visible(False) # Turn off y-axis labels
    plt.show()
    return
def snpy_histogram_plotter(path, raw=False, save_plot=False, save_loc='../default/', ignore_type=[], param_bins=[None, None, None, None, None]):
    # Get reviewed fits
    with open(gen.get_constants()['reviewed_fits_txt'], 'r') as f:
        reviewed_good_fits, reviewed_okay_fits, reviewed_bad_fits = [], [], []
        for container in [reviewed_good_fits, reviewed_okay_fits, reviewed_bad_fits]:
            for i in f.readline().split(', ')[1:]:
                container.append(i[1:-1])
            container[-1] = container[-1][:-1]

    # Pull data
    objs = gen.dict_handler(path=path, choice='unpack')
    mu, st, Tmax, EBVhost, host_mass = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
    for obj in objs:
        # Clean data
        if not raw:
            objname = obj
            if objname[:2] == 'SN':
                objname = objname[2:]
            if ('good' in ignore_type) and (objname in reviewed_good_fits):
                continue
            if ('okay' in ignore_type) and (objname in reviewed_okay_fits):
                continue
            if ('bad' in ignore_type) and (objname in reviewed_bad_fits):
                continue

        mu = np.append(mu, float(objs[obj]['mu']))
        st = np.append(st, float(objs[obj]['st']))
        Tmax = np.append(Tmax, float(objs[obj]['Tmax']))
        EBVhost = np.append(EBVhost, float(objs[obj]['EBVhost']))
        host_mass = np.append(host_mass, float(objs[obj]['host_mass']))

    # Plot
    fig, ax = plt.subplots(1, 5, figsize=(16, 4), layout='constrained')
    params = [mu, st, Tmax, EBVhost, host_mass]
    param_names = ['mu', 'st', 'Tmax', 'EBVhost', 'host_mass']
    param_bins = [45, 45, 45, 45, 45]
    for i in range(len(params)):
        ax[i].hist(params[i], bins=param_bins[i])
        if i != 0:
            ax[i].get_yaxis().set_visible(False)
        ax[i].set_xlabel(param_names[i])

    plt.suptitle("Parameters for '" + path.split('/')[-1].split('_')[0].upper()
                 + "' data\n Number of Transients: " + str(len(objs)), fontsize=20)
    if save_plot:
        print('Saved figure to... ', save_loc+path.split('/')[-1].split('_')[0]+'_hist.png')
        plt.savefig(save_loc+path.split('/')[-1].split('_')[0]+'_hist.png')
    plt.show()

    return
def salt3_histogram_plotter(path, raw=False, save_plot=False, save_loc='../default/', ignore_type=[], param_bins=[None, None, None, None, None]):
    # Pull data
    objs = gen.dict_handler(path=path, choice='unpack')
    t0, x0, x1, c, mu, host_mass = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
    for obj in objs:
        t0 = np.append(t0, float(objs[obj]['t0']))
        x0 = np.append(x0, float(objs[obj]['x0']))
        x1 = np.append(x1, float(objs[obj]['x1']))
        c = np.append(c, float(objs[obj]['c']))
        mu = np.append(mu, float(objs[obj]['mu']))
        host_mass = np.append(host_mass, float(objs[obj]['host_mass']))

    # Plot
    fig, ax = plt.subplots(1, 6, figsize=(25, 5), layout='constrained')
    params = [t0, x0, x1, c, mu, host_mass]
    param_names = ['t0', 'x0', 'x1', 'c', 'mu', 'host_mass']
    param_bins = [45, 45, 45, 45, 45, 45]
    for i in range(len(params)):
        ax[i].hist(params[i], bins=param_bins[i])
        if i != 0:
            ax[i].get_yaxis().set_visible(False)
        ax[i].set_xlabel(param_names[i])
        ax[i].set_xlabel(param_names[i])
        ax[i].set_xlabel(param_names[i])

    plt.suptitle("Parameters for '" + path.split('/')[-1].split('_')[0].upper()
                 + "' data\n Number of Transients: " + str(len(objs)), fontsize=20)
    if save_plot:
        print('Saved figure to... ', save_loc+path.split('/')[-1].split('_')[0]+'_hist.png')
        plt.savefig(save_loc+path.split('/')[-1].split('_')[0]+'_hist.png')
    plt.show()

    return