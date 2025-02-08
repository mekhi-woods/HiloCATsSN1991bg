import os
import requests
import numpy as np
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM

def get_apikeys(apikeys_loc: str = 'txts/api_keys.txt') -> dict:
    """
    :param apikeys_loc: Location of API keys file.
    :return: dictionary with API keys.
    """
    # Check file exsists
    if not os.path.isfile(apikeys_loc):
        raise FileNotFoundError('[!!!] API keys file not found!')

    # Read api_keys.txt
    APIKEY = {}
    with open(apikeys_loc, 'r') as f:
        # Check if file empty
        temp = f.readlines()
        if len(temp) > 0:
            for line in temp:
                line = line[:-1].split(', ')
                if len(line) != 2:
                    continue
                APIKEY.update({line[0]: line[1]})

    # If empty write user input to api_keys.txt
    if len(APIKEY) == 0:
        with open(apikeys_loc, 'w') as f:
            print('[!!!] API keys file empty! Please enter the following API keywords...')
            f.write(f"tns_bot_id, {input('tns_bot_id: ')}\n")
            f.write(f"tns_bot_name, {input('tns_bot_name: ')}\n")
            f.write(f"tns_bot_api_key, {input('tns_bot_api_key: ')}\n")
            f.write(f"atlas_key, {input('atlas_key [Optional, Press enter to skip...]: ')}\n")
    return APIKEY
def get_constants(constant_loc='constants.txt'):
    CONSTANTS = {}
    with open(constant_loc, 'r') as f:
        temp = f.readlines()
        for line in temp:
            line = line[:-1].split(', ')
            if len(line) != 2:
                continue
            CONSTANTS.update({line[0]: line[1]})
    return CONSTANTS
def current_cosmo(H0=70, O_m=0.3):
    return FlatLambdaCDM(H0, O_m)
def default_open(path: str, table_mode: bool = False):
    """
    Open file that has the typical format of this project's 'txt' files.
    :param path: str; location of file.
    :param table_mode: bool; whether or not to return the data as an astropy table.
    :return: (list, np.array[str]) | astropy.table.Table.
    """
    data = np.genfromtxt(path, dtype='str', delimiter=', ')
    hdr, data = data[0, :], data[1:, :]
    for i in range(len(hdr)):
        if hdr[i] in ['objname', 'origin', 'algo']: continue
        data[:, i] = data[:, i].astype(float)
    if table_mode:
        var_table = Table()
        hdr = hdr.tolist()
        for h in hdr:
            try: var_table[h] = data[:, hdr.index(h)].astype(float)
            except ValueError: var_table[h] = data[:, hdr.index(h)]
        return var_table
    else:
        return hdr.tolist(), data
def get_twomass():
    if not os.path.exists('twomass++_velocity_LH11.npy'):
        raise FileNotFoundError('[!!!] TwoMass velocities file "twomass++_velocity_LH11.npy" not found! '
                                'Please download from the follow...\n'
                                "https://drive.usercontent.google.com/download?id=1DGcWQPgmI2ZoHJm_zCqyscogmWwY7lQu&"
                                "export=download&authuser=0")
    return