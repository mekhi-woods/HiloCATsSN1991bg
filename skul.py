import warnings
with warnings.catch_warnings():
    import scripts.general as gen
    from astropy.cosmology import FlatLambdaCDM
    import matplotlib.pyplot as plt
    import numpy as np
    import sys
    import os

CURRENT_COSMO = FlatLambdaCDM(70, 0.3)  # Hubble Constant, Omega-Matter
CONSTANTS = gen.get_constants()

def norm_Ia_objs_pull(loc='../txts/normIaparams.txt', quiet=False):
    # Check quiet
    if quiet:
        sys.stdout = open(os.devnull, 'w')

    normIaObjs = {}
    with open(loc) as f:
        hdr = []
        for n in f.readline().split(' '):
            if n != '':
                hdr.append(n)
        hdr[-1] = hdr[-1][:-1]
        print(hdr)

        data = f.readlines()
        i, t = 0, len(data[:-1])
        for line in data[:-1]:
            i += 1
            print('---------------------------------------------------------------------------------------------------')
            print('[ '+str(i)+' / '+str(t)+' ]')
            n_line = line.split(' ')
            n_line[-1] = n_line[-1][:-1] # Remove '\n' from last element
            n_line[0] = n_line[0][:-4] # Remove the .txt from ID
            ra, dec = float(n_line[1]), float(n_line[2]) # Save RA/DEC

            objname, z = gen.TNS_objname_z(ra, dec) # Save Objname

            normIaObjs.update({objname: {
                'ra': ra, 'dec': dec, 'z': z, 'mu': n_line[4], 'mu_err': 0.00,
                't0': n_line[5], 'x0': n_line[6], 'x1': n_line[7], 'c': n_line[8],
                't0err': n_line[9], 'x0err': n_line[10], 'x1err': n_line[11], 'cerr': n_line[12] # 'original_name': n_line[0],
            }})

    # Restore print statements
    sys.stdout = sys.__stdout__

    return normIaObjs
def sort_resid_z(objs, x1lim=1):
    # objnames, all_resid, all_resid_err, all_z = [], [], [], []
    # for obj in objs:
    #     # Normalize Data
    #     if str(objs[obj]['mu']) == 'nan' or str(objs[obj]['z']) == 'nan':
    #         continue
    #     elif str(objs[obj]['mu']) == 'None' or str(objs[obj]['z']) == 'None':
    #         continue
    #     # Last adjustments
    #     elif float(objs[obj]['x1']) < -x1lim or float(objs[obj]['x1']) > x1lim:
    #         continue
    #     else:
    #         z = float(objs[obj]['z'])
    #         resid = float(objs[obj]['mu']) - CURRENT_COSMO.distmod(float(objs[obj]['z'])).value
    #         resid_err = float(objs[obj]['mu_err'])
    #
    #         all_z.append(z)
    #         all_resid.append(resid)
    #         all_resid_err.append(resid_err)
    #         objnames.append(obj)
    # return objnames, all_resid, all_resid_err, all_z




    return
def plot_z_v_mu(objnames, all_resid, all_resid_err, all_z, title, x1lim='', labels=False):
    fig, axs = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'width_ratios': [10, 1]}, constrained_layout=True)

    # Plotting
    axs[0].errorbar(all_z, all_resid, yerr=all_resid_err, fmt='o')
    axs[1].hist(all_resid, bins=40, orientation="horizontal")
    if labels:
        for i in range(len(objnames)):
            axs[0].text(all_z[i], all_resid[i], objnames[i], fontsize=8)

    # Formatting
    ylimiter = np.max(np.abs(all_resid))+0.5
    axs[0].set_ylim(-ylimiter, ylimiter); axs[1].set_ylim(-ylimiter, ylimiter)
    fig.suptitle(title +  # Figure Title
                 '-' + str(x1lim) + ' < x1 < ' + str(x1lim) + ' | Median: ' + str(round(np.median(all_resid), 2)) + ' | '
                 'Scatter: ' + str(round(np.std(all_resid), 2)) + ' | # of pts: ' + str(len(all_resid)), size='medium')
    axs[0].set(xlabel='Redshift', ylabel='Hubble Residuals')  # Sub-plot Labels
    axs[1].get_yaxis().set_visible(False) # Turn off y-axis labels
    plt.show()
    return
def plot_mass_v_mu(objnames, all_resid, all_z, title, x1lim='', labels=False):
    fig, axs = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'width_ratios': [10, 1]}, constrained_layout=True)

    # Plotting
    axs[0].scatter(all_z, all_resid)
    axs[1].hist(all_resid, bins=40, orientation="horizontal")
    if labels:
        for i in range(len(objnames)):
            axs[0].text(all_z[i], all_resid[i], objnames[i], fontsize=8)

    # Formatting
    ylimiter = np.max(np.abs(all_resid))+0.5
    axs[0].set_ylim(-ylimiter, ylimiter); axs[1].set_ylim(-ylimiter, ylimiter)
    fig.suptitle(title +  # Figure Title
                 '-' + str(x1lim) + ' < x1 < ' + str(x1lim) + ' | Median: ' + str(round(np.median(all_resid), 2)) + ' | '
                 'Scatter: ' + str(round(np.std(all_resid), 2)) + ' | # of pts: ' + str(len(all_resid)), size='medium')
    axs[0].set(xlabel='Redshift', ylabel='Hubble Residuals')  # Sub-plot Labels
    axs[1].get_yaxis().set_visible(False) # Turn off y-axis labels
    plt.show()
    return

if __name__ == '__main__':
    # params = {'x1lim': 0.5}
    #
    # # Plot Normal SNe Ia
    # # txtNormPath = '../default/sneNormObjs.txt'
    #
    # # sneNormObjs = norm_Ia_objs_pull(quiet=True)
    # # objnames, all_resid, all_resid_err, all_z = sort_resid_z(sneNormObjs, **params)
    # # gen.dict_handler(choice='pack', data_dict=sneNormObjs, path=txtNormPath)
    #
    # # gen.host_mass(txtNormPath, save_loc='../default/', keep_data=False, update_saved=True, use_mass_key=True)
    #
    # # plot_z_v_mu(objnames, all_resid, all_resid_err, all_z, title='Hubble Residuals vs. Redshift of Normal SNe Ia\n', labels=True, **params)
    #
    #
    # # Plot 91bg-like SNe Ia
    # # txt91bgPath = '../default/sne91bgObjs.txt'
    # #
    # # sn91bgIaObjs = gen.dict_handler(choice='unpack', path='../saved/salt/atlas/salt_atlas_saved.txt')
    # # print(len(sn91bgIaObjs))
    # # objnames, all_resid, all_resid_err, all_z = sort_resid_z(sn91bgIaObjs, **params)
    # #
    # # gen.dict_handler(choice='pack', data_dict=sn91bgIaObjs, path=txt91bgPath)
    # gen.host_mass(txt91bgPath, save_loc='../default/', keep_data=False, update_saved=True, use_mass_key=True)
    #
    # plot_z_v_mu(objnames, all_resid, all_resid_err, all_z, title='Hubble Residuals vs. Redshift of 91bg-like SNe Ia\n', labels=True, **params)

    # gen.mass_step_calc(txtNormPath)

    objs = gen.dict_handler(choice='unpack', path=CONSTANTS['atlas_saved_loc'] + 'atlas_saved.txt')
