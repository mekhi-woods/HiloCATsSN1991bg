import os
import json
import glob
import astropy
import requests
import numpy as np
import time as systime
import matplotlib.pyplot as plt
from astropy.table import Table
from collections import OrderedDict
from astropy.time import Time as astrotime

class sneObj:
    # Construction Functions ----------------------------------------------------------------------------------------- #
    def __init__(self, source: str, path: str):
        source = source.lower()
        if source == 'empty' or len(path) == 0:
            print('[+++] Creating empty SNe object...')
            pass
        elif source == 'class' and os.path.exists(path):
            print('[+++] Creating SNe object from class file...')
            self.load_class(path)
        elif source in ['atlas-91bg','atlas-norm', 'csp-91bg', 'csp-norm', 'ztf-91bg', 'ztf-norm'] and os.path.exists(path):
            print(f"[+++] Creating '{source}' SNe object using '{path}'...")

            # Load and clean data
            var_tbl = self.make_objTbl(source, path)
            var_tbl = self.clean_objTbl(var_tbl)
            self.get_details(np.average(var_tbl['ra']), np.average(var_tbl['dec']))

            # Set arrays
            self.mjd = np.array(var_tbl['mjd'])
            self.mag = np.array(var_tbl['mag'])
            self.dmag = np.array(var_tbl['dmag'])
            self.flux = np.array(var_tbl['flux'])
            self.dflux = np.array(var_tbl['dflux'])
            self.filter = np.array(var_tbl['filter'])
            self.zp = np.array(var_tbl['zp'])
            self.params = {} # Initalized for later
            self.covariance = [] # Initalized for later

            # Extra details
            self.origin = source
            self.path = f"classes/{source}/{source}_{self.objname}_class.txt"
            if source == 'csp' and self.z == None: # CSP seems to sometimes fail to find a 'z' with TNS
                try: self.z = round(np.average(var_tbl['z']), 4)
                except: pass

            # Save class file
            self.save_class()
        return
    def __str__(self):
        return (f"{self.objname}, {self.origin} @ {self.path}\n"
                f"{str(self.coords)} | ({self.z}, {str(round(self.z_cmb, 2))})")
    def make_objTbl(self, source: str, path: str) -> astropy.table.Table:
        # Load data and make intial table
        if source == 'atlas-91bg':
            # Open data
            with open(path, 'r') as f: hdr = f.readline()[1:-1].split(',')
            data = np.genfromtxt(path, dtype='str', delimiter=',')

            # Make table
            var_table = Table()
            for h in hdr:
                try: var_table[h] = data[:, hdr.index(h)].astype(float)
                except ValueError: var_table[h] = data[:, hdr.index(h)]

            # Add parity with CSP & ZTF
            var_table.remove_columns(['expname', 'flux', 'dflux', 'snr', 'snr', 'x', 'y', 'peakval', 'skyval', 'peakfit',
                                      'dpeak', 'skyfit', 'chin', 'major', 'minor', 'snrdet', 'snrlimit', 'apfit',
                                      'mag5sig', 'texp'])
            for h_old, h_new in zip(['dm', 'uJy', 'duJy'],
                                    ['dmag', 'flux', 'dflux']):
                var_table[h_old].name = h_new
        elif source == 'atlas-norm':
            # Load data
            with open(path, 'r') as f: hdr = f.readline()[1:-1].split(',')
            data = np.genfromtxt(path, dtype='str', delimiter=',', skip_header=1)

            # Make table
            var_table = Table()
            for h in hdr:
                try: var_table[h] = data[:, hdr.index(h)].astype(float)
                except ValueError: var_table[h] = data[:, hdr.index(h)]

            # Add parity with CSP & ZTF
            var_table.remove_columns(['err', 'chi/N', 'x', 'y', 'maj', 'min', 'phi', 'apfit', 'mag5sig', 'Sky', 'Obs'])
            var_table['zp'] = np.full(len(var_table), np.nan)
            for h_old, h_new in zip(['JD', 'm', 'dm', 'uJy', 'duJy', 'F', 'RA', 'Dec'],
                                    ['mjd', 'mag', 'dmag', 'flux', 'dflux', 'filter', 'ra', 'dec']):
                var_table[h_old].name = h_new
        elif source == 'csp-91bg' or source == 'csp-norm':
            var_table = Table(names=['filter', 'zp', 'z', 'ra', 'dec', 'mjd', 'mag', 'dmag', 'flux', 'dflux'],
                              dtype=[str, float, float, float, float, float, float, float, float, float])
            with open(path, 'r') as f:
                csp_objname, csp_z, csp_ra, csp_dec = f.readline()[2:-1].split(' ')
                for l in f.readlines():
                    l = l.split(' ')
                    # Filter line
                    if len(l) == 2:
                        csp_filter = str(l[1][:-1])
                        csp_zp = float(get_constant()['csp_zpts_'+csp_filter])
                    else:
                        csp_mjd, csp_mag, csp_dmag = float(l[-3])+53000, float(l[-2]), float(l[-1])
                        csp_flux = 10 ** ((csp_mag - csp_zp) / -2.5)
                        csp_dflux = np.abs(csp_flux) * np.log(10) * ((1 / 2.5) * csp_dmag)

                        var_table.add_row([csp_filter, csp_zp, csp_z, csp_ra, csp_dec, csp_mjd,
                                           csp_mag, csp_dmag, csp_flux, csp_dflux])
        elif source == 'ztf-91bg' or source == 'ztf-norm':
            # Load hdr & data
            with open(path, 'r') as f:
                for i in range(3): f.readline()
                ztf_ra = float(f.readline().split(' ')[-2])
                ztf_dec = float(f.readline().split(' ')[-2])
            data = np.genfromtxt(path, delimiter=' ', dtype=str, skip_header=54)
            hdr, data = list(data[0]), data[1:]
            for i in range(len(hdr)): hdr[i] = hdr[i][:-1]

            # Make table
            var_table = Table()
            for h in hdr:
                try:
                    var_table[h] = data[:, hdr.index(h)].astype(float)
                except ValueError:
                    var_table[h] = data[:, hdr.index(h)]

            # Add RA & DEC
            var_table['ra'] = np.full(len(var_table), ztf_ra)
            var_table['dec'] = np.full(len(var_table), ztf_dec)

            # Fix time, JD to MJD
            var_table['jd'] = np.array(var_table['jd']).astype(float) - 2400000.5
            var_table['forcediffimflux'][var_table['forcediffimflux'] == 'null'] = 'nan'
            var_table['forcediffimfluxunc'][var_table['forcediffimfluxunc'] == 'null'] = 'nan'
            var_table['mag'] = (-2.5 * np.log10(np.array(var_table['forcediffimflux']).astype(float))) + var_table['zpdiff']
            var_table['dmag'] = np.abs(-1.08573620476 * (np.array(var_table['forcediffimfluxunc']).astype(float)
                                                         / np.array(var_table['forcediffimflux']).astype(float)))

            # Add parity with CSP & ATLAS
            var_table.remove_columns(['index', 'field', 'ccdid', 'qid', 'pid', 'infobitssci', 'sciinpseeing',
                                      'scibckgnd', 'scisigpix', 'zpmaginpsci', 'zpmaginpsciunc', 'zpmaginpscirms',
                                      'clrcoeff', 'clrcoeffunc', 'ncalmatches', 'exptime', 'adpctdif1', 'adpctdif2',
                                      'diffmaglim', 'programid', 'rfid', 'forcediffimsnr', 'forcediffimchisq',
                                      'forcediffimfluxap', 'forcediffimfluxuncap', 'forcediffimsnrap', 'aperturecorr',
                                      'dnearestrefsrc', 'nearestrefmag', 'nearestrefmagunc', 'nearestrefchi',
                                      'nearestrefsharp', 'refjdstart', 'refjdend', 'procstatu'])
            for h_old, h_new in zip(['zpdiff', 'jd', 'forcediffimflux', 'forcediffimfluxunc'],
                                    ['zp', 'mjd', 'flux', 'dflux']):
                var_table[h_old].name = h_new
        else:
            raise ValueError(f'[!!!] Unknown source, {source}! Must be [atlas-91bg/atlas-norm/csp/ztf]...')
        return var_table
    def clean_objTbl(self, tbl: astropy.table.Table) -> astropy.table.Table:
        # Remove '>' from ATLAS mags
        for i, n_mag in enumerate(tbl['mag']):
            if str(n_mag)[0] == '>': tbl['mag'][i] = float(n_mag[1:])

        # Remove nulls
        for col in ['mag', 'dmag', 'flux', 'dflux']:
            tbl = tbl[~np.isnan(np.array(tbl[col]).astype(float))]
        return tbl
    def get_details(self, ra: float, dec: float):
        success, r = False, 2
        print(f"[+++] Querying TNS @ ({ra}, {dec}), with radius={r}...")
        while not success:
            try:
                details = get_TNSDetails(ra, dec, radius=r)
                self.objname = details['objname']
                self.z = details['redshift']
                self.z_cmb = np.nan
                self.coords = [details['radeg'], details['decdeg']]
                self.discdate = float(astrotime(details['discoverydate'], format='iso').jd) - 2400000.5  # JD to MJD
                success = True
            except:
                print(f"[~~~] Warning: TNS found no targets @ ({ra}, {dec}), with radius={r}... Raising to r={r+1}")
                r += 1
            if r == 10: raise RuntimeError('[!!!!!] TNS could not find transient!')
        return
    def save_class(self):
        print(f"[+++] Saving {self.objname} class to {self.path}...")
        with open(self.path, 'w') as f:
            f.write(f"{self.origin},{self.objname},{self.coords[0]},{self.coords[1]},{self.z},{self.z_cmb},{self.discdate}\n")
            f.write('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
            for p in self.params:
                f.write(f"{p},{self.params[p]['value']},{self.params[p]['err']}\n")
            f.write('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
            for line in self.covariance:
                p_line = str(line[0])
                for i in range(1, len(line)):
                    p_line += ',' + str(line[i])
                f.write(p_line + '\n')
            f.write('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
            for col, name in zip([self.zp, self.filter, self.mjd, self.flux, self.dflux, self.mag, self.dmag],
                                 ['zp', 'filter', 'mjd', 'flux', 'dflux', 'mag', 'dmag']):
                f.write(f"{name}")
                for n_col in col:
                    f.write(f",{n_col}")
                f.write('\n')
    def load_class(self, path: str):
        print(f"[+++] Opening class from {path}...")
        with open(path, 'r') as f:
            # Read header
            details = f.readline().split(',')
            self.path = path
            self.origin = details[0]
            self.objname = details[1]
            self.coords = [float(details[2]), float(details[3])]
            self.z = float(details[4])
            self.z_cmb = float(details[5])
            self.discdate = float(details[6][:-1])
            f.readline()  # Skip break line

            # Read params
            line = f.readline()
            while '+++' not in line:
                line = line.split(',')
                self.params.update({line[0]: {'value': float(line[1]), 'err': float(line[2])}})
                line = f.readline()

            # Read covariances
            if 'salt' in path:
                self.covariance = np.array([])
                line = f.readline()
                while '+++' not in line:
                    line = line.split(', ')
                    line[-1] = line[-1][:-1]
                    self.covariance = np.append(self.covariance, np.array(line).astype(float))
                    line = f.readline()
                self.covariance = self.covariance.reshape(4, 4)
            f.readline()  # Skip break line

            # Read arrays
            self.zp = np.array(f.readline().split(',')[1:])
            self.zp[-1] = self.zp[-1][:-1]
            self.filter = np.array(f.readline().split(',')[1:])
            self.filter[-1] = self.filter[-1][:-1]
            self.mjd = np.array(f.readline().split(',')[1:])
            self.mjd[-1] = self.mjd[-1][:-1]
            self.flux = np.array(f.readline().split(',')[1:])
            self.flux[-1] = self.flux[-1][:-1]
            self.dflux = np.array(f.readline().split(',')[1:])
            self.dflux[-1] = self.dflux[-1][:-1]
            self.mag = np.array(f.readline().split(',')[1:])
            self.mag[-1] = self.mag[-1][:-1]
            self.dmag = np.array(f.readline().split(',')[1:])
            self.dmag[-1] = self.dmag[-1][:-1]
        return

    # Display Functions ---------------------------------------------------------------------------------------------- #
    def plot(self, y_type: str = 'mag', save_loc: str = ''):
        print(f"[+++] Plotting LC of '{self.objname}' ({y_type})...")
        filter_dict = {'u': 'teal', 'g': 'green', 'r': 'red', 'i': 'indigo', 'B': 'blue',
                       'V0': 'violet', 'V1': 'purple', 'V': 'red', 'Y': 'goldenrod', 'Hdw': 'tomato', 'H': 'salmon',
                       'J': 'aquamarine', 'Jrc2': 'cadetblue', 'Jdw': 'turquoise', 'Ydw': 'olive',
                       'c': 'cyan', 'o': 'orange', 'ZTF_g': 'green', 'ZTF_r': 'red', 'ZTF_i': 'indigo'}

        # Select y-axis element
        y_axis = self.mag if y_type == 'mag' else self.flux
        y_axis_err = self.dmag if y_type == 'mag' else self.dflux

        # Plot
        fig, axs = plt.subplots(1, 1, figsize=(12, 6), constrained_layout=True)
        for f in np.unique(self.filter):
            axs.errorbar(self.mjd[self.filter == f].astype(float), y_axis[self.filter == f].astype(float),
                         yerr=y_axis_err[self.filter == f].astype(float),
                         color = filter_dict[f], fmt='o', ms=3)

        # Lines
        axs.axvline(self.discdate, color='black', linestyle='--', label='Discovery Date')
        try: axs.axvline(self.params['Tmax']['value'], color='maroon', linestyle='--', label='Peak Brightness')
        except: pass

        # Formatting
        if y_type == 'mag': axs.invert_yaxis()
        axs.legend()
        axs.set_xlabel('MJD', size=16)
        axs.set_ylabel('Magnitude', size=16) if y_type == 'mag' else axs.set_ylabel('Flux (uJy)', size=16)
        plt.suptitle(f"Lightcurve of '{self.objname}' ({y_type})", size=16)
        if len(save_loc) != 0:
            print('[+++] Saving to '+save_loc)
            plt.savefig(save_loc)
        plt.show()
        return

    # Fitting Functions ---------------------------------------------------------------------------------------------- #

def get_constant(cosntant_loc='txts/constants.txt'):
    CONSTANTS = {}
    with open(cosntant_loc, 'r') as f:
        temp = f.readlines()
        for line in temp:
            line = line[:-1].split(', ')
            if len(line) != 2:
                continue
            CONSTANTS.update({line[0]: line[1]})
    return CONSTANTS
def get_APIkeys(apikeys_loc='api_keys.txt'):
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
def get_TNSDetails(ra, dec, radius=2):
    APIKEY = get_APIkeys()
    tns_bot_id, tns_bot_name, tns_bot_api_key = APIKEY['tns_bot_id'], APIKEY['tns_bot_name'], APIKEY['tns_bot_api_key']
    tns_marker = (
        f'tns_marker{{"tns_id": "{int(tns_bot_id)}",'
        f'"type": "bot", "name": "{tns_bot_name}"}}'
    )
    headers = {"User-Agent": tns_marker}
    search_obj = [
        ("ra", str(ra)),
        ("dec", str(dec)),
        ("radius", str(radius)),
        ("units", "arcsec"),
        ("objname", ""),
        ("objname_exact_match", 0),
        ("internal_name", ""),
        ("internal_name_exact_match ", 0),
        ("objid", ""),
        ("public_timestamp", ""),
    ]
    search_data = {"api_key": tns_bot_api_key, "data": json.dumps(OrderedDict(search_obj))}
    response = requests.post("https://www.wis-tns.org/api/get/search", headers=headers, data=search_data)
    response = json.loads(response.text)
    transients =  response["data"]
    get_obj = [
        ("objname", transients[0]["objname"]),
        ("objid", transients[0]["objid"]),
        ("photometry", "0"),
        ("spectra", "0"),
    ]
    get_data = {"api_key": tns_bot_api_key, "data": json.dumps(OrderedDict(get_obj))}
    response = requests.post("https://www.wis-tns.org/api/get/object", headers=headers, data=get_data)
    response = json.loads(response.text)
    details = response["data"]
    return details
def fit(data_loc: str):
    paths = glob.glob(data_loc)
    dataset = paths[0].split('/')[-2].lower()

    # Verify proper dataset
    if dataset not in ['atlas-91bg','atlas-norm', 'csp-91bg', 'csp-norm', 'ztf-91bg', 'ztf-norm']:
        raise ValueError(f"[!!!] Dataset, '{dataset}', not recognized!")

    # Select fitting mode
    if len(paths) == 1:
        print(f"[+++] Fitting file '{data_loc}'...")  # Indivisual fit
        sn = sneObj(dataset, paths[0])  # Initiate

        return sn
    elif len(paths) > 1:
        print(f"[+++] Fitting data in '{data_loc}'...")  # Batch fit
        sne = []
        for i, path in enumerate(paths):
            print(f'[{i + 1} / {len(paths)}] ================================================================')
            sne.append(sneObj(dataset, path))  # Initiate
        return sne
    else:
        print('[!!!] Invalid file/data path!')
        return



    # all_atlas_norm = glob.glob('data/ATLAS-norm/*.txt')
    # for i, path in enumerate(all_atlas_norm):
    #     print(f'[{i + 1} / {len(all_atlas_norm)}] ====================================================================')
    #     test = sneObj('atlas-norm', path)
    # all_atlas_91bg = glob.glob('data/ATLAS-91bg/*.txt')
    # for i, path in enumerate(all_atlas_91bg):
    #     print(f'[{i + 1} / {len(all_atlas_91bg)}] ====================================================================')
    #     test = sneObj('atlas-91bg', path)
    #
    # all_csp_norm = glob.glob('data/CSP-norm/*.txt')
    # for i, path in enumerate(all_csp_norm):
    #     print(f'[{i + 1} / {len(all_csp_norm)}] ====================================================================')
    #     test = sneObj('csp', path)
    #
    # all_csp_91bg = glob.glob('data/CSP-91bg/*.txt')
    # for i, path in enumerate(all_csp_91bg):
    #     print(f'[{i + 1} / {len(all_csp_91bg)}] ====================================================================')
    #     test = sneObj('csp', path)
    #
    # all_ztf_91bg = glob.glob('data/ZTF-91bg/*.txt')
    # for i, path in enumerate(all_ztf_91bg):
    #     print(f'[{i + 1} / {len(all_ztf_91bg)}] ====================================================================')
    #     test = sneObj('ztf', path)
    return


if __name__ == '__main__':
    start = systime.time()  # Runtime tracker
    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')
