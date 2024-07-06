import os
import numpy as np
import scripts.general as gen
import requests
import json
import matplotlib.pyplot as plt

def ztf_wget(submit=False):
    # ra, dec, jds, jde = 186.07860833333334, 10.446222222222222, 2460350, 2460450
    # email = 'mekhidw@hawaii.edu'  # email you subscribed with.
    # userpass = 'wxdk286'  # password that was issued to you.
    # cmd = f"wget --http-user=ztffps --http-passwd=dontgocrazy! -O ../tests/log.txt \"https://ztfweb.ipac.caltech.edu/cgi-bin/requestForcedPhotometry.cgi?ra={dec}&dec={ra}&jdstart={jds}&jdend={jde}&email={'mekhidw@hawaii.edu'}&userpass={'wxdk286'}\""
    # print(cmd)
    # os.system(cmd)

    atlas_txt = gen.dict_unpacker('../snpy/atlas/atlas_saved.txt')

    print('Requesting the following from ZTF...')
    for tar in atlas_txt:
        print('---------------------------------------')
        ra = atlas_txt[tar]['ra']
        dec = atlas_txt[tar]['dec']
        jds = float(atlas_txt[tar]['MJDs'])
        jde = float(atlas_txt[tar]['MJDe'])

        if submit:
            cmd = f"wget --http-user=ztffps --http-passwd=dontgocrazy! -O ../tests/log.txt \"https://ztfweb.ipac.caltech.edu/cgi-bin/requestForcedPhotometry.cgi?ra={dec}&dec={ra}&jdstart={jds}&jdend={jde}&email={'mekhidw@hawaii.edu'}&userpass={'wxdk286'}\""
            # print(cmd)
            results = os.system(cmd)
            print(results)

            print('RA: '+str(ra)+', DEC: '+str(dec)+', '+str(jds)+'-'+str(jde)+', Status: ')
        break



    return

def ztf_processing():
    header_num = 54
    path = '../data/ZTF/old/batchfp_req0002064117_lc.txt'
    data = np.genfromtxt(path, delimiter=' ', skip_header=header_num, dtype=str)

    # Grab header
    hdr = []
    with open(path, 'r') as f:
        for i in range(header_num):
            hdr.append(f.readline())

    # Create Object
    # flux(forcediffimflux), flux_unc(forcediffimfluxunc)
    obj = {'ra': hdr[3].split(' ')[5], 'dec': hdr[4].split(' ')[5],
           'time': data[:, 22].astype(float), 'flux': data[:, 24].astype(float), 'flux_unc': data[:, 25].astype(float)}

    # Clean data
    negs = np.where(obj['flux'] > 0)[0]
    extremes = np.where(obj['flux_unc'] < 150)[0]


    fixes = np.unique(np.concatenate((negs, extremes)))

    for cat in ['flux', 'flux_unc', 'time']:
        obj[cat] = obj[cat][negs]

    sigma = 1
    plt.errorbar(obj['time'], obj['flux'], yerr=obj['flux_unc']*sigma, fmt='o')
    # plt.gca().invert_yaxis()
    plt.show()


    return