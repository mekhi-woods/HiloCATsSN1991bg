import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import matplotlib.pyplot as plt
    import snpy
    import glob

files = glob.glob('../snpy/ztf/ascii/*_snpy.txt')
for i in range(len(files)):
    print('['+str(i)+']\n------------------------------------------------------------------------------------------')

    # n_s = snpy.get_sn('../snpy/ztf/ascii/2020acoo_snpy.txt')
    n_s = snpy.get_sn(files[i])
    n_s.choose_model('EBV_model2', stype='st')

    # Remove empty filters -- fix for 'ValueError: attempt to get argmin of an empty sequence'
    for filter in list(n_s.data.keys()):
        if len(n_s.data[filter].magnitude) == 0:
            del n_s.data[filter]

    # Fitting
    run = True
    while run:
        try:
            n_s.fit(bands=None, dokcorr=True, k_stretch=True, reset_kcorrs=True)
            print('Success')
            run = False

            n_s.plot()
            plt.show()

            print('Results:\n',
                  '\t mu =', n_s.parameters['DM'], '+/-', n_s.errors['DM'], '\n',
                  '\t st =', n_s.parameters['st'], '+/-', n_s.errors['st'], '\n',
                  '\t Tmax =', n_s.parameters['Tmax'], '+/-', n_s.errors['Tmax'], '\n',
                  '\t EBVhost =', n_s.parameters['EBVhost'], '+/-', n_s.errors['EBVhost'], '\n')
        except Exception as error:
            if 'All weights for filter' and 'are zero.' in str(error):
                print('Weights for filter', str(error)[23:24], 'are zero. Removing...')
                del n_s.data[str(error)[23:24]]
            elif 'Error:  to solve for EBVhost, you need to fit more than one filter' in str(error):
                print('Only one filter, cannot fit data.')
                run = False
            else:
                print(error)
                run = False






