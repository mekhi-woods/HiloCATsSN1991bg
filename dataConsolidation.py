import os
import glob
import utils  # Import of utils.py
import datetime
import numpy as np
import time as systime
from astropy.table import Table

CURRENTDATE = datetime.datetime.now()

def make_param_file(sne: list, save_loc: str):
    with open(save_loc, 'w') as f:
        print(f"# Created by M.D. Woods -- {CURRENTDATE} -- NUM TARGETS: {len(sne)}", file=f)
        print(f"# WARNING: For files with SNPY & SALT, the 'stretch' and 'color' are NOT identical.", file=f)
        for sn in sne:
            print(sn)


    return

if __name__ == '__main__':
    start = systime.time()  # Runtime tracker
    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')
