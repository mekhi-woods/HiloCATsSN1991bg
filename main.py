# M. D. Woods
# 02/08/2025
import time as systime
import utils              # Utility functions for various common actions accross project; Here to be initiated for later.
import queryCSP           # Retrieves "Photometric Data release 3" tarball from Krisciunas et al. 2017 hosted on
                          # https://csp.obs.carnegiescience.edu/data, then seperates data using target file.
import queryATLAS         # Retrieves data from ATLAS server and stores it in data\ -- requires 'atlas_key' in 'txts/api_keys.txt'.
import queryZTF           # Retrieves data from ZTF server sends out request to be send to an email.
import fitter             # Uses data in data\ to create SNooPy & SALT fits.
import plotter            # Plots data found in files in results\.
import dataConsolidation  # Consolidates list of sneObjs into one TXT parameter file.

def run_queryCSP():
    queryCSP.download(save_loc='data/')  # location of CSP-norm & CSP-91bg
    queryCSP.seperate_data(tar_list='txts/target_files/tarlist_CSP_91bg.csv',
                           data_loc='data/',
                           save_loc='data/CSP-91bg')  # Seperate 1991bg-like out of all CSP-I data
    queryCSP.seperate_data(tar_list='txts/target_files/tarlist_CSP_norm.csv',
                           data_loc='data/',
                           save_loc='data/CSP-norm')  # Seperate Normals out of all CSP-I data
    return
def run_queryATLAS():
    queryATLAS.download(tar_list = "txts/target_files/tarlist_ATLAS_91bg.csv", save_loc = 'data/ATLAS-91bg/')
    queryATLAS.download(tar_list = "txts/target_files/tarlist_ATLAS_norm.csv", save_loc = 'data/ATLAS-norm/')
    return
def run_queryZTF():
    queryZTF.download(tns_targets_path = "txts/target_files/tarlist_ZTF_91bg.csv", save_loc = 'data/ZTF-91bg/')
    queryZTF.download(tns_targets_path = "txts/target_files/tarlist_ZTF_norm.csv", save_loc = 'data/ZTF-norm/')
    return
def run_fitter(algo: str = 'snpy', rewrite: bool = True):
    """
    :param algo: Algorithm to use for fitting
    :param rewrite: Whether to rewrite the fit data or not
    :return: None
    """
    # Batch
    sne_csp_91bg = fitter.fit(data_loc = 'data/CSP-91bg/*.txt', algo = algo, rewrite = rewrite)
    # sne_csp_norm = fitter.fit(data_loc = 'data/CSP-norm/*.txt', algo = algo, rewrite = rewrite)
    # sne_atlas_91bg = fitter.fit(data_loc = 'data/ATLAS-91bg/*.txt', algo = algo, rewrite = rewrite)
    # sne_atlas_norm = fitter.fit(data_loc = 'data/ATLAS-norm/*.txt', algo = algo, rewrite = rewrite)
    # sne_ztf_91bg = fitter.fit(data_loc = 'data/ZTF-91bg/*.txt', algo = algo, rewrite = rewrite)
    # sne_ztf_norm = fitter.fit(data_loc = 'data/ZTF-norm/*.txt', algo = algo, rewrite = rewrite)  # Data not downloaded

    ## Combine SNe into one param file
    dataConsolidation.make_param_file(sne_csp_91bg, f'results/csp_91bg_{algo}.txt')
    # dataConsolidation.make_param_file(sne_csp_91bg+sne_csp_norm+sne_atlas_91bg+sne_atlas_norm+sne_ztf_91bg)

    # Indivisual
    # print("================================================================")
    # fitter.fit(data_loc = 'data/CSP-91bg/SN2007ax_snpy.txt', algo = algo, rewrite = rewrite)
    # print("================================================================")
    # fitter.fit(data_loc = 'data/CSP-norm/SN2005be_snpy.txt', algo = algo, rewrite = rewrite)
    # print("================================================================")
    # fitter.fit(data_loc = 'data/ATLAS-91bg/ATLAS2020abmg.txt', algo = algo, rewrite = rewrite)
    # print("================================================================")
    # fitter.fit(data_loc = 'data/ATLAS-norm/ATLAS2024zjo.txt', algo = algo, rewrite = rewrite)
    # print("================================================================")
    # fitter.fit(data_loc = 'data/ZTF-91bg/forcedphotometry_req00381085_lc.txt', algo = algo, rewrite = rewrite)
    # print("================================================================")
    # fitter.fit(data_loc = 'data/ZTF-norm/*.txt', algo = algo)  # None downloaded yet

    # Testing
    # fitter.fit(data_loc = 'data/ATLAS-norm/ATLAS2024msu.txt', algo = algo, rewrite = False)
    # sne = fitter.fit(data_loc = 'data/CSP-91bg/SN2006bd_snpy.txt', algo = algo, rewrite = True)

    return
def run_plotter():
    final_dir = 'plots/'

    pm_norm_salt_cut = 'results/old/aaronDo_salt2_params_cut.txt'
    pm_norm_salt_uncut = 'results/old/aaronDo_salt2_params.txt'
    pm_norm_snpy_cut = 'results/old/dr3_params.txt'  # Needs to be cut?, no mu
    pm_norm_snpy_uncut = 'results/old/dr3_params.txt'
    pm_norm_merged_cut = 'results/old/aaronDo_salt2_params_cut.txt'  # only contains salt fitted
    pm_norm_merged_uncut = 'results/old/aaronDo_salt2_params.txt'  # only contains salt fitted

    pm_91bg_salt_cut = 'results/old/combiend__salt_params_cut.txt'
    pm_91bg_salt_uncut = 'results/old/combiend__salt_params.txt'
    pm_91bg_snpy_cut = 'results/old/combiend__snpy_params_cut.txt'
    pm_91bg_snpy_uncut = 'results/old/combiend__snpy_params.txt'
    pm_91bg_merged_cut = 'results/old/merged_params_cut.txt'
    pm_91bg_merged_uncut = 'results/old/merged_params.txt'

    pm_redNorms = 'results/old/redNormSNe_salt.txt'
    pm_dust = 'results/old/global_dust_params.txt'

    plotter.resid_v_mass(path_91bg=pm_91bg_merged_cut,
                          path_norm=pm_norm_merged_cut,
                          save_loc=final_dir+'resid_v_mass.png',
                          label = False)
    plotter.mu_v_z(path_91bg=pm_91bg_merged_cut,
                    path_norm=pm_norm_merged_cut,
                    save_loc=final_dir+'mu_v_z.png',
                    label = False)

    ## SALT3 Plots
    plotter.alpha_beta(path_91bg=pm_91bg_salt_cut,
                        path_norm=pm_norm_salt_cut,
                        save_loc=final_dir+'alpha_beta.png')

    ## Dust Plots
    plotter.abs_mag_v_dust(path_91bg=pm_91bg_salt_cut, path_red_norm=pm_redNorms, path_dust=pm_dust,
                            save_loc=final_dir+'absMag_v_dust.png')
    plotter.dust_hist(path_91bg=pm_91bg_salt_cut, path_red_norm=pm_redNorms, path_dust=pm_dust,
                       save_loc=final_dir+'dust_params.png')
    plotter.resid_v_mass_dust(path_91bg=pm_91bg_merged_cut, path_norm=pm_norm_merged_cut, path_dust=pm_dust,
                               save_loc=final_dir+'dust_resid_v_mass.png')

    ## Paramater Histograms
    plotter.param_hist(snpy_91bg_path=pm_91bg_snpy_uncut, snpy_norm_path=pm_norm_snpy_uncut,
                       salt_91bg_path=pm_91bg_salt_uncut, salt_norm_path=pm_norm_salt_uncut,
                       save_loc=final_dir + 'param_hist_uncut.png', line_type='median')
    plotter.param_hist(snpy_91bg_path=pm_91bg_snpy_cut, snpy_norm_path=pm_norm_snpy_cut,
                       salt_91bg_path=pm_91bg_salt_cut, salt_norm_path=pm_norm_salt_cut,
                       save_loc=final_dir + 'param_hist_cut.png', line_type='median')

    ## Brout+Scolnic 2021 style Color v. Scatter
    plotter.color_v_scatter(path_snpy_91bg=pm_91bg_snpy_cut, path_snpy_norm='results/norm_snpy_params_cut.txt',
                            path_salt_91bg=pm_91bg_salt_cut, path_salt_norm=pm_norm_salt_cut,
                            bin_nums=[[10, 10], [10, 10]], bin_bounds=[[-0.3, 0.6], [-0.3, 0.6]], label=True,
                            save_loc=final_dir + 'color_v_scatter.png')

    ## Brout+Scolnic 2021 style Dust v. Scatter
    plotter.dust_v_scatter(pm_91bg_merged_cut, pm_dust, bin_num=60, bin_bounds = [0, 7], label=True,
                           save_loc=final_dir + 'dust_v_scatter.png')
    return

if __name__ == '__main__':
    start = systime.time()  # Runtime tracker
    # utils.get_twomass()  # Checks for twomass velocities data (large file)
    # run_queryCSP()
    # run_queryATLAS()
    # run_queryZTF()
    run_fitter(algo = 'snpy', rewrite = False)
    # run_fitter(algo = 'salt', rewrite = True)
    # run_plotter()
    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')
