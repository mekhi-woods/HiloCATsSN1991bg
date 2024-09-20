import os
import sys
import glob
import shutil
import time as systime
import numpy as np

# from astro_ghost.ghostHelperFunctions import getTransientHosts
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from astropy.time import Time as astrotime

from astroquery.sdss import SDSS
from astroquery.mast import Catalogs
import warnings

def current_cosmo(H0=70, O_m=0.3):
    return FlatLambdaCDM(H0, O_m)
def get_constants(cosntant_loc='txts/constants.txt'):
    CONSTANTS = {}
    with open(cosntant_loc, 'r') as f:
        temp = f.readlines()
        for line in temp:
            line = line[:-1].split(', ')
            if len(line) != 2:
                continue
            CONSTANTS.update({line[0]: line[1]})
    return CONSTANTS
def get_APIkeys(apikeys_loc='txts/APIkeys.txt'):
    if not os.path.isfile(apikeys_loc):
        raise FileNotFoundError('[!!!] API keys file not found!')
    APIKEY = {}
    with open(apikeys_loc, 'r') as f:
        temp = f.readlines()
        for line in temp:
            line = line[:-1].split(', ')
            if len(line) != 2:
                continue
            APIKEY.update({line[0]: line[1]})
    return APIKEY
def TNS_details(obj_ra, obj_dec, attempts=10, radius=2, use_keys=True):
    import tnsAPI
    APIkey = get_APIkeys()
    print('-----------------------------------------------------------------------------------------------------------')
    print('Searching TNS for ['+str(obj_ra)+', '+str(obj_dec)+']...')
    tns_bot_id, tns_bot_name, tns_bot_api_key = APIkey['tns_bot_id'], APIkey['tns_bot_name'], APIkey['tns_bot_api_key']
    tns_key = dict_handler(choice='unpack', path=get_constants()['tns_key_txt'])
    obj_name, obj_z, obj_discdate = '', 0.00, ''

    # Check key
    print('Checking TNS key...')
    if use_keys and (str(obj_ra) in tns_key):
        obj_name, obj_z, obj_discdate = tns_key[str(obj_ra)]['objname'], tns_key[str(obj_ra)]['z'], tns_key[str(obj_ra)]['discoverydate']
    else:
        # Query TNS key
        print('Not in TNS key, querying TNS...')
        try:
            # Code abridged from David's code
            headers = tnsAPI.build_tns_header(tns_bot_id, tns_bot_name)
            tns_api_url = f"https://www.wis-tns.org/api/get"
            search_tns_url = tnsAPI.build_tns_url(tns_api_url, mode="search")
            get_tns_url = tnsAPI.build_tns_url(tns_api_url, mode="get")
            search_data = tnsAPI.build_tns_search_query_data(tns_bot_api_key, obj_ra, obj_dec, radius=radius)
            transients = tnsAPI.rate_limit_query_tns(search_data, headers, search_tns_url)
            get_data = tnsAPI.build_tns_get_query_data(tns_bot_api_key, transients[0])
            details = tnsAPI.rate_limit_query_tns(get_data, headers, get_tns_url)
            obj_name, obj_z = details['objname'], details['redshift']
            obj_discdate = float(astrotime(details['discoverydate'], format='iso').jd) - 2400000.5 # JD to MJD
            if use_keys:
                with open(get_constants()['tns_key_txt'], 'a') as f:
                    f.write(str(obj_ra) + ', ' + str(obj_dec) + ', ' + obj_name + ', ' + str(obj_z) +
                            ', ' + str(obj_discdate) + '\n') # Save to key
        except Exception as error:
            print('Error:', error)
            if attempts <= 0:
                raise RuntimeError('TNS completly timed out.')
            print(f'{attempts} attempts remaining')
            systime.sleep(5)
            obj_name, obj_z, obj_discdate = TNS_details(obj_ra, obj_dec, attempts=attempts-1, radius=radius+1, use_keys=use_keys)
    return obj_name, obj_z, obj_discdate
def dict_handler(data_dict={}, choice='', path='default/dict.txt', delimiter=', '):
    if choice == 'unpack':
        print('[+++] Unpacking objects from '+path+'...')
        with open(path, 'r') as f:
            hdr = f.readline()[:-1].split(delimiter)
        data = np.genfromtxt(path, delimiter=delimiter, dtype=str, skip_header=1)
        if len(data) == 0:
            return {}
        temp_objs = {}

        if len(np.shape(data)) == 1:
            iter_period = len(data[0])
        else:
            iter_period = len(data[:, 0])

        for i in range(iter_period):
            obj = data[:, 0][i]
            temp_objs.update({obj: {}})
            for j in range(len(hdr)):
                temp_objs[obj].update({hdr[j]: data[i, j]})
        return temp_objs
    elif choice == 'pack':
        print('[+++] Packing objects into text file...')
        if len(data_dict) == 0:
            print('[!!!] Attempted to pack an empty dictionary!')
            return
        catagories = list(data_dict[list(data_dict.keys())[0]].keys())
        with open(path, 'w') as f:
            f.write('objname')
            for category in catagories:
                f.write(delimiter+category)
            f.write('\n')
            for objname in data_dict:
                f.write(objname)
                for category in catagories:
                    f.write(delimiter+str(data_dict[objname][category]))
                f.write('\n')
        print('[+++] Files packed to', path)
        return
    elif choice == 'arrays':
        print('[+++] Unpacking objects from '+path+' into arraies...')
        with open(path, 'r') as f:
            hdr = f.readline()[:-1].split(delimiter)
        data = np.genfromtxt(path, delimiter=delimiter, dtype=str, skip_header=1)
        if len(data) == 0:
            return {}
        temp_arrs = {}
        for i in range(len(hdr)):
            temp_arrs.update({hdr[i]: data[:, i]})
        return temp_arrs
    else:
        print("[!!!] Invalid packing option! ['arrays'/'pack'/'unpack']")
        return


