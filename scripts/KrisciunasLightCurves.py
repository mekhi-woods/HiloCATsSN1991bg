import matplotlib.pyplot as plt
import numpy as np
import time
import os

FILTER_WHEEL = ['u', 'g', 'r', 'i', 'B', 'V0']
CUR_TIME = str(int(time.time()))

def KrLC(save=False, saveLoc='save\\'):
    if save:
        workingDir = saveLoc + '\\' + CUR_TIME
        os.mkdir(workingDir)

    KrisciunasPath = "targetLists\\91bglike_justnames.txt"
    KrisciunasNames = np.genfromtxt(KrisciunasPath, dtype=str, delimiter=', ')

    allCPSPhot = "data\\CSPdata\\SN_photo.dat"
    allCPSPhotData = np.genfromtxt(allCPSPhot, dtype='str')

    names = allCPSPhotData[:, 0]
    filters = allCPSPhotData[:, 1]
    time = allCPSPhotData[:, 2]
    light = allCPSPhotData[:, 3]
    err = allCPSPhotData[:, 4]

    sigma = 1
    for tar in KrisciunasNames:
        for n in range(len(FILTER_WHEEL)):
            # output_names = names[(names == tar) & (filters == FILTER_WHEEL[n])]
            output_light = light[(names == tar) & (filters == FILTER_WHEEL[n])].astype('float64')
            output_time = time[(names == tar) & (filters == FILTER_WHEEL[n])].astype('float64') + 53000
            output_err = err[(names == tar) & (filters == FILTER_WHEEL[n])].astype('float64')
            plt.errorbar(output_time, output_light, yerr=output_err*sigma, fmt='o', label=FILTER_WHEEL[n])

        plt.title(tar); plt.xlabel('Time [MJD]'); plt.ylabel('Intensity [mag]')
        plt.gca().invert_yaxis()
        plt.legend()
        if save:
            name = str(tar)+'.png'
            plt.savefig(workingDir+'\\'+name)
            print("Saved '"+name+"' @ "+workingDir)
        plt.show()


    # print(CUR_TIME)

    # Alternative Format
    # sigma = 1
    # for tar in KrisciunasNames:
    #     for n in range(len(FILTER_WHEEL)):
    #         # output_names = names[(names == tar) & (filters == FILTER_WHEEL[n])]
    #         output_light = light[(names == tar) & (filters == FILTER_WHEEL[n])].astype('float64')
    #         output_time = time[(names == tar) & (filters == FILTER_WHEEL[n])].astype('float64')
    #         output_time_mod = output_time - np.min(output_time)
    #         output_err = err[(names == tar) & (filters == FILTER_WHEEL[n])].astype('float64')
    #         plt.errorbar(output_time_mod, output_light, yerr=output_err * sigma, fmt='o', label=FILTER_WHEEL[n])
    #
    #     plt.title(tar)
    #     # noinspection PyUnboundLocalVariable
    #     plt.xlabel('MJD - '+str(np.min(output_time+53000)))
    #     plt.ylabel('Intensity [mag]')
    #     plt.gca().invert_yaxis()
    #     plt.legend()
    #     if save:
    #         plt.savefig(saveLoc+str(tar)+'.png')
    #     plt.show()

    return