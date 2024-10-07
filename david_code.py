# import astropy.table as at
# import matplotlib.pyplot as plt
# import numpy as np
# import sncosmo
# from astropy.stats import sigma_clipped_stats
# # from scipy.stats import binned_statistic
# from astropy.cosmology import FlatLambdaCDM
#
# if __name__ == '__main__':
#     # data = at.Table.read('output/combiend__snpy_params.txt' ,format='ascii')
#     # data = at.Table.read('output/combiend__snpy_params_91bgkver.txt' ,format='ascii')
#     data = at.Table.read('output/combiend__snpy_params_david.txt' ,format='ascii')
#
#

#


import astropy.table as at
import matplotlib.pyplot as plt
import numpy as np
import sncosmo
from astropy.stats import sigma_clipped_stats
# from scipy.stats import binned_statistic
from astropy.cosmology import FlatLambdaCDM

if __name__ == '__main__':
    data = at.Table.read('output/combiend__snpy_params.txt' ,format='ascii')
    # data = at.Table.read('output/combiend__snpy_params_91bgkver.txt' ,format='ascii')
    # data = at.Table.read('output/combiend__snpy_params_david.txt' ,format='ascii')

    data = data[(data['EBVhost'] > -0.2) & (data['EBVhost'] < 0.3)]
    data = data[data['st_err'] < 0.1]
    data = data[data['EBVhost_err'] < 0.1]
    data = data[data['Tmax_err'] < 1]
    data = data[data['mu_err'] < 0.2]
    data = data[(data['st'] < 1.0)]
    data = data[(data['z'] > 0.015)]
    resid = data['mu'] - FlatLambdaCDM(70, 0.3).distmod(data['z_cmb']).value
    print(sigma_clipped_stats(resid)[2])
