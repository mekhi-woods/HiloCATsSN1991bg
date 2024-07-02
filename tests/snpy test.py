import warnings



import matplotlib.pyplot as plt

n_s = snpy.get_sn('../data/CSPdata/SN2005bo_snpy.txt')
n_s.fit(bands=None, dokcorr=True, k_stretch=False, reset_kcorrs=True)
n_s.plot()
plt.show()
# print(n_s.__dict__)

# print('Results:',
#       'mu =', n_s.parameters['DM'], '+/-', n_s.errors['DM'], '\n',
#       '\t st =', n_s.parameters['st'], '+/-', n_s.errors['st'], '\n',
#       '\t Tmax =', n_s.parameters['Tmax'], '+/-', n_s.errors['Tmax'], '\n',
#       '\t EBVhost =', n_s.parameters['EBVhost'], '+/-', n_s.errors['EBVhost'], '\n')