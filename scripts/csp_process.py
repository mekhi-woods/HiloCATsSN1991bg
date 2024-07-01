import numpy as np

def burns_fitting(skip_problems=False, use_saved=False, snpy_plots=True):
    print('[+++] Fitting CSP data with SNooPy...')
    # Get Chris Burns Data
    objnames = np.genfromtxt('/content/HiloCATsSN1991bg/targetLists/91bglike_justnames.txt', dtype=str, delimiter=', ')

    # Get CSP paths of Chris Burns Data
    objpaths = []
    for name in objnames:
        if os.path.exists('/content/HiloCATsSN1991bg/data/CSPdata/SN'+name+'_snpy.txt'):
            objpaths.append('/content/HiloCATsSN1991bg/data/CSPdata/SN'+name+'_snpy.txt')
        else:
            print(name, 'not found...')

    # Fitting
    objs = {}
    err_i = 0
    fit_args = {'skip_problems': skip_problems, 'use_saved': use_saved, 'snpy_plots': snpy_plots}
    print('Fitting arguments: ', fit_args)
    for n in range(len(objpaths)):
        tracker = '['+str(n+1)+'/'+str(len(objpaths))+']' # Purely cosmetic
        objname = objpaths[n][39:-9]
        print(tracker, objname)
        temp_dict = snpy_fit(objpaths[n], objname, save_loc=SNPY_BURNS, plot_save_loc=SNPY_BURNS_PLOTS, **fit_args)

        if type(temp_dict) == dict:
            print('\tResults: \n',
                  '\t\tmu =', temp_dict['mu'], '+/-', temp_dict['mu_err'],'\n',
                  '\t\tst =', temp_dict['st'], '+/-', temp_dict['st_err'],'\n',
                  '\t\tTmax =', temp_dict['Tmax'], '+/-', temp_dict['Tmax_err'],'\n',
                  '\t\tEBVhost =', temp_dict['EBVhost'],  '+/-', temp_dict['EBVhost_err'])
            print('\tMJD min:', temp_dict['MJDs'], '+/- 100 MJD', '| MJD max:', temp_dict['MJDe'], '+/- 100 MJD\n')
            objs.update({objname: temp_dict})
        else:
            err_i += 1
            print('[!!!]\t', temp_dict, '\n')
    print('Finshed! Successfully fit', len(objpaths)-err_i, '/', len(objpaths), 'SNe from ATLAS! ['+str(round(((len(objpaths)-err_i) / len(objpaths))*100, 2))+'%]')

    # Announce problem children
    print('Problem children: ', handle_problem_children(state='READ'))

    # Save Data
    if len(objs) > 0:
        dict_packer(objs, BURNS_SAVE_TXT, delimiter=', ') # Save data from fitting
    return

def burns_plotting(choice):
    if 'reg_hist' in choice:
        print('[+++] Ploting Histogram of SNooPy Fitting Parameters...')
        data = np.genfromtxt(BURNS_SAVE_TXT, dtype=str, delimiter=', ', skip_header=1)
        objnames, st, Tmax, EBVHost = data[:, 0], data[:, 3].astype(float), data[:, 4].astype(float), data[:, 5].astype(float)
        reg_params = [st, Tmax, EBVHost]
        bins_reg = [15, 20, 15]
        plot_title = 'CSP SNooPy Parameters'

        fig, ax = plt.subplots(1, 3, figsize=(15, 5))
        titles = ['st', plot_title+'\nTmax', 'EBVhost']
        for i in range(len(titles)):
            ax[i].hist(reg_params[i], bins=bins_reg[i])
            ax[i].set_title(titles[i])
        plt.show()

    if 'res_hist' in choice:
        print('[+++] Ploting Histogram of SNooPy-Chris Burns Parameters Residuals...')
        data = np.genfromtxt(BURNS_SAVE_TXT, dtype=str, delimiter=', ', skip_header=1)
        objnames, st, Tmax, EBVHost = data[:, 0], data[:, 3].astype(float), data[:, 4].astype(float), data[:, 5].astype(float)

        # Take the difference between CSP and DR3 file
        st_res, Tmax_res, EBVHost_res = [], [], []
        dr3 = read_DR3()
        for n in range(len(objnames)):
            if objnames[n] in dr3:
                st_res.append(st[n] - dr3[objnames[n]]['st'])
                Tmax_res.append(Tmax[n] - dr3[objnames[n]]['Tmax'])
                EBVHost_res.append(EBVHost[n] - dr3[objnames[n]]['EBVHost'])\

        # Correct for MJD -- 53000
        for i in range(len(Tmax_res)):
            Tmax_res[i] = Tmax_res[i] + 53000

        # Plot
        res_params = [st_res, Tmax_res, EBVHost_res]
        bins_res = [30, 50, 20]
        xlims = [[-0.1, 0.1], [-3, 3], [-0.15, 0.15]]
        plot_title = 'SNooPy-Chris Burns Parameters'
        fig, ax = plt.subplots(1, 3, figsize=(15, 5))
        titles = ['st', plot_title+'\nTmax', 'EBVhost']
        for i in range(len(titles)):
            ax[i].hist(res_params[i], bins=bins_res[i])
            ax[i].set_title(titles[i])
            ax[i].set_xlim(xlims[i])

        plt.show()

    if 'zvmu' in choice:
        print('[+++] Ploting redshift [z] vs distance mod [mu] of SNooPy Parameters...')
        data = np.genfromtxt(BURNS_SAVE_TXT, dtype=str, delimiter=', ', skip_header=1)
        objnames, mu, z = data[:, 0], data[:, 1].astype(float), data[:, 2].astype(float)

        plt.figure(figsize=(8, 4))
        for n in range(len(objnames)):
            plt.loglog(z[n], mu[n], label=objnames[n], marker='o')
        plt.title('CSP Redshift vs. Distance Mod\n SNe # = '+str(len(objnames)))
        plt.xlabel('Redshift')
        plt.ylabel('Distance Mod')
        plt.legend()
        plt.show()

    return
