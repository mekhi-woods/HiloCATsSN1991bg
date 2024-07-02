import os
import snpy
import numpy as np
import matplotlib.pyplot as plt
from scripts import tns_redshifts

""" GENERAL """
def TNS_details(ra, dec):
    tns_bot_id, tns_bot_name, tns_bot_api_key = '73181', 'YSE_Bot1', '0d771345fa6b876a5bb99cd5042ab8b5ae91fc67'

    # Code from David
    headers = tns_redshifts.build_tns_header(tns_bot_id, tns_bot_name)
    tns_api_url = f"https://www.wis-tns.org/api/get"

    # get the API URLs
    search_tns_url = tns_redshifts.build_tns_url(tns_api_url, mode="search")
    get_tns_url = tns_redshifts.build_tns_url(tns_api_url, mode="get")

    search_data = tns_redshifts.build_tns_search_query_data(tns_bot_api_key, ra, dec)
    transients = tns_redshifts.rate_limit_query_tns(search_data, headers, search_tns_url)

    get_data = tns_redshifts.build_tns_get_query_data(tns_bot_api_key, transients[0])
    transient_detail = tns_redshifts.rate_limit_query_tns(get_data, headers, get_tns_url)

    return transient_detail
def snpy_fit(path, objname, save_loc, plot_save_loc,
             skip_problems=False, use_saved=False, snpy_plots=True, save_plots=True):
    problem_children = handle_problem_children(state='READ') # Open problem children

    if skip_problems and (objname in problem_children):
        return None
    else:
        try:
            if use_saved and os.path.exists(save_loc+objname+'_EBV_model2.snpy'):
                n_s = snpy.get_sn(save_loc+objname+'_EBV_model2.snpy')
            else:
                n_s = snpy.get_sn(path)
                n_s.choose_model('EBV_model2', stype='st')
                n_s.set_restbands() # Auto pick appropriate rest-bands

                # Sort out empty filters & get start and end time
                mjds, mjde = [], []
                filter_wheel = []
                for filter in list(n_s.data.keys()):
                    if len(n_s.data[filter].MJD) <= 3:
                        print('\t', objname, 'has too few points in', filter, 'filter')
                        continue
                    mjds.append(min(n_s.data[filter].MJD))
                    mjde.append(max(n_s.data[filter].MJD))
                    filter_wheel.append(filter)
                print('\t', filter_wheel)

                n_s.fit(bands=filter_wheel, dokcorr=True, k_stretch=False, reset_kcorrs=True, **{'mangle':1,'calibration':0})
                n_s.save(save_loc+objname+'_EBV_model2.snpy')

            if snpy_plots:
                n_s.plot(outfile=plot_save_loc+objname+'_snpyplots.png')
                plt.show()
        except Exception as error:
            problem_children = np.append(problem_children, objname)
            handle_problem_children(state='WRITE', problem_c=problem_children) # Commit problem children
            return error

    plt.close()
    return {'ra': n_s.ra, 'dec': n_s.decl, 'z': n_s.z, 'MJDs': min(mjds), 'MJDe': max(mjde),
            'mu': n_s.parameters['DM'], 'st': n_s.parameters['st'], 'Tmax': n_s.parameters['Tmax'], 'EBVhost': n_s.parameters['EBVhost'],
            'mu_err': n_s.errors['DM'], 'st_err': n_s.errors['st'], 'Tmax_err': n_s.errors['Tmax'], 'EBVhost_err': n_s.errors['EBVhost']}
def snpy_fit_indv(objname):
    path = '/content/HiloCATsSN1991bg/snpy/atlas/'+objname+'_EBV_model2.snpy'

    n_s = snpy.get_sn(path)
    n_s.choose_model('EBV_model2', stype='st')
    n_s.set_restbands() # Auto pick appropriate rest-bands

    # Sort out empty filters & get start and end time
    mjds, mjde = [], []
    filter_wheel = []
    for filter in list(n_s.data.keys()):
        if len(n_s.data[filter].MJD) <= 3:
            print('\t', objname, 'has too few points in', filter, 'filter')
            continue
        mjds.append(min(n_s.data[filter].MJD))
        mjde.append(max(n_s.data[filter].MJD))
        filter_wheel.append(filter)

    n_s.fit(bands=filter_wheel, dokcorr=True, k_stretch=False, reset_kcorrs=True, **{'mangle':1,'calibration':0})
    plt.show()

    print('Results:',
        'mu =', n_s.parameters['DM'], '+/-', n_s.errors['DM'],'\n',
        '\t st =', n_s.parameters['st'], '+/-', n_s.errors['st'],'\n',
        '\t Tmax =', n_s.parameters['Tmax'], '+/-', n_s.errors['Tmax'],'\n',
        '\t EBVhost =', n_s.parameters['EBVhost'],  '+/-', n_s.errors['EBVhost'],'\n',
        '\t MJD min:', min(mjds), '| MJD max:', max(mjde))

    return
def write_ASCII(objs, save_loc, quiet=True):
    print('[+++] Saving data to ASCII files for SNooPy...')
    for obj in objs:
        with open(save_loc+obj+'_snpy.txt', 'w') as f:
            # Line 1 -- Objname, Helio-Z, RA, Dec (Ex. SN1981D 0.005871 50.65992 -37.23272)
            f.write(str(obj)+' '+str(objs[obj]['z'])+' '+str(objs[obj]['ra'])+' '+str(objs[obj]['dec'])+'\n')

            # 'o'/'ATri'-filter photometry block -- Date (JD/MJD), mag, err (674.8593 12.94 0.11)
            f.write('filter ATri\n')
            for i in range(len(objs[obj]['time_o'])):
                f.write(str(objs[obj]['time_o'][i])+'\t'+str(objs[obj]['mag_o'][i])+'\t'+str(objs[obj]['dmag_o'][i])+'\n')

            # # 'c'/'ATgr'-filter photometry block
            f.write('filter ATgr\n')
            for i in range(len(objs[obj]['time_c'])):
                f.write(str(objs[obj]['time_c'][i])+'\t'+str(objs[obj]['mag_c'][i])+'\t'+str(objs[obj]['dmag_c'][i])+'\n')
    return
def read_DR3(loc='../txts/DR3_fits.dat'):
    data = np.genfromtxt(loc, dtype=str, skip_header=1)
    dr3 = {}
    for n in range(len(data[:, 0])):
        dr3.update({data[:, 0][n]: {'st': float(data[:, 1][n]), 'e_st': float(data[:, 2][n]), 'z': float(data[:, 3][n]),
                           'Tmax': float(data[:, 5][n]), 'e_Tmax': float(data[:, 6][n]),
                           'EBVHost': float(data[:, 25][n]), 'e_EBVHost': float(data[:, 26][n])}})
    return dr3
def dict_unpacker(path, delimiter=', '):
    with open(path, 'r') as f:
        hdr = f.readline()[:-1].split(delimiter)

    data = np.genfromtxt(path, delimiter=delimiter, dtype=str, skip_header=1)
    if len(data) == 0:
        return {}
    temp_objs = {}
    for i in range(len(data[:, 0])):
        obj = data[:, 0][i]
        temp_objs.update({obj: {}})
        for j in range(len(hdr)):
            temp_objs[obj].update({hdr[j]: data[i, j]})
    return temp_objs
def dict_packer(data_dict, save_loc, delimiter=', '):
    catagories = list(data_dict[list(data_dict.keys())[0]].keys())
    with open(save_loc, 'w') as f:
        f.write('objname')
        for category in catagories:
            f.write(delimiter+category)
        f.write('\n')
        for objname in data_dict:
            f.write(objname)
            for category in catagories:
                f.write(delimiter+str(data_dict[objname][category]))
            f.write('\n')
    return
def handle_problem_children(state, problem_c=None):
    path = '../problem_children.txt'

    if state == 'READ':
        # Read problem children
        problem_c = np.genfromtxt(path, dtype=str)
        return problem_c
    elif state == 'WRITE':
        # Write problem children
        problem_c = np.unique(problem_c)
        with open(path, 'w') as f:
            for c in problem_c:
                f.write(c+'\n')
        return None
    else:
        raise Exception("Invalid state: '"+state+"' [READ/WRITE]")
def save_to_zip(zip_loc, save_loc):
    print('Saving zipped files to...', save_loc)
    files = glob.glob(zip_loc+'*')
    with ZipFile(save_loc, 'w') as zip:
        for n_file in files:
            zip.write(n_file)
    return
