import time as systime
import queryATLAS  # Retrieves data from ATLAS server and stores it in \data\
import fitter      # Uses data in \data\ to create SNooPy & SALT fits
import plotter     # Plots data found in files in \results\

def run_queryATLAS():
    queryATLAS.download(tns_targets_path = "txts/tns_targets_normIa.csv", save_loc = 'data/ATLAS-norm/')
    return
def run_fitter():
    fitter.fit('data/ATLAS-norm/*.txt')
    # fitter.fit('data/ATLAS-norm/ATLAS2024afze.txt')
    return
def run_plotter():
    final_dir = 'plots/'

    pm_norm_salt_cut = 'results/aaronDo_salt2_params_cut.txt'
    pm_norm_salt_uncut = 'results/aaronDo_salt2_params.txt'
    pm_norm_snpy_cut = 'results/dr3_params.txt'  # Needs to be cut?, no mu
    pm_norm_snpy_uncut = 'results/dr3_params.txt'
    pm_norm_merged_cut = 'results/aaronDo_salt2_params_cut.txt'  # only contains salt fitted
    pm_norm_merged_uncut = 'results/aaronDo_salt2_params.txt'  # only contains salt fitted

    pm_91bg_salt_cut = 'results/combiend__salt_params_cut.txt'
    pm_91bg_salt_uncut = 'results/combiend__salt_params.txt'
    pm_91bg_snpy_cut = 'results/combiend__snpy_params_cut.txt'
    pm_91bg_snpy_uncut = 'results/combiend__snpy_params.txt'
    pm_91bg_merged_cut = 'results/merged_params_cut.txt'
    pm_91bg_merged_uncut = 'results/merged_params.txt'

    pm_redNorms = 'results/redNormSNe_salt.txt'
    pm_dust = 'results/global_dust_params.txt'

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
    # run_queryATLAS()
    # run_fitter()
    # run_plotter()
    print('|---------------------------|\n Run-time: ', round(systime.time() - start, 4), 'seconds\n|---------------------------|')
