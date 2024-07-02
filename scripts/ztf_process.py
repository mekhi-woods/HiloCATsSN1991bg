import os
import numpy as np
import scripts.general as gen
import requests
import json

def ztf_collection(submit=False):
    limit = 10
    ra, dec, jds, jde = [], [], [], []
    atlas_txt = gen.dict_unpacker('../snpy/atlas/atlas_saved.txt')

    for tar in atlas_txt:
        ra.append(float(atlas_txt[tar]['ra']))
        dec.append(float(atlas_txt[tar]['dec']))
        jds.append(float(atlas_txt[tar]['MJDs']))
        jde.append(float(atlas_txt[tar]['MJDe']))

    ra_json, dec_json = json.dumps(ra), json.dumps(dec)
    # jdstart, jdend = json.dumps(min(jds)), json.dumps(max(jde))
    jdstart, jdend = 2458216.1234, 2458450.0253

    if submit:
        payload = {'ra': ra_json, 'dec': dec_json, 'jdstart': jdstart, 'jdend': jdend, 'email': 'mekhidw@hawaii.edu', 'userpass': 'wxdk286'}
        r = requests.post('https://ztfweb.ipac.caltech.edu/cgi-bin/batchfp.py/submit', auth=('ztffps', 'dontgocrazy!'), data=payload)
        print("Status_code=", r.status_code)
        print(r.text)



    return

def ztf_processing():
    header_num = 54
    path = '../data/ZTF/batchfp_req0002064117_lc.txt'
    data = np.genfromtxt(path, delimiter=' ', skip_header=header_num, dtype=str)

    hdr = []
    with open(path, 'r') as f:
        for i in range(header_num):
            hdr.append(f.readline())

    obj = {'ra': hdr[3].split(' ')[5], 'dec': hdr[4].split(' ')[5]}

    print(data[:, 22]) #JD
    print(data[:, 10]) # mag (zpmaginpsci)
    print(data[:, 11]) # magerr (zpmaginpsciunc)




    return