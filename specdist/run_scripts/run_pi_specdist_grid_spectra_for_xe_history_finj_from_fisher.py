import specdist as sd
import numpy as np

finj_from_fisher = 'no'
add_edges = 'yes'
photon_injection_case = 'lyc_reio'

# load photon injection and cosmotherm modules
ct = sd.cosmotherm()
cosmo = sd.cosmo()
X_dm = sd.dm_particle()



sd_lib = sd.specdist_ct_spectra_lib()
sd.load_ct_spectra_lib(photon_injection_case,sd_lib)


firas = sd.firas()
edges = sd.edges()
a_dict = {}
a_dict['firas'] = firas
a_dict['edges'] = edges
a_dict['add_edges'] = 'yes'

# set relevant parameter values
ct.ct_Drho_rho_dec = 1e-300
# ct.ct_h = 0.70
# ct.ct_Omega_b = 0.0457
# ct.ct_Omega_m = 0.30
ct.ct_emission_absorption_mode = 0
ct.ct_npts = 5000
ct.ct_zend = 1.1e6
ct.ct_zstart = 2e6

ct.ct_include_pi = 1
ct.ct_evolve_Xe = 1
ct.save_Xe = 'yes'
ct.ct_lyc = 1
ct.ct_reionisation_model = 1

sd.set_dm_params_to_CT_pi_params(X_dm,ct)
sd.set_cosmo_to_CT_cosmo_params(cosmo,ct)

xi_array = np.logspace(np.log10(1e-8),np.log10(1.e7),64)
#xi_array = np.logspace(np.log10(1e-8),np.log10(1.e7),64)
xi_array = [1e-6]
#Gamma_X = 1e-6 #np.logspace(-9,-6,4)
#zi_array = [1.1e6,9.9e5]
for Gamma_X in [np.logspace(-12,-17,10)[0]]:
    str_gamma = str("%.3e"%Gamma_X)

    ct.ct_Gamma_dec = Gamma_X
    ct.ct_x_dec = 1e7

    if finj_from_fisher == 'yes':
        Gamma_values = [ct.ct_Gamma_dec]
        xinj_array = [ct.ct_x_dec]
        f_dm_fisher = sd.pi_run_fisher_constraints(Gamma_values,xinj_array,sd_lib,**a_dict)
        print(f_dm_fisher)
        ct.ct_photon_injection_f_dec = f_dm_fisher['curves'][0]['finj']
        print('finj_fisher = %e'%ct.ct_photon_injection_f_dec)
    else:
        ct.ct_photon_injection_f_dec = 0.
        print('not using finj')



    args = {}
    p_name = 'photon injection x_dec'
    #p_name = 'zend'
    p_array = xi_array
    #p_array = zi_array
    args['param_values_array'] = p_array
    args['param_name'] = p_name
    args['save_spectra'] = 'yes'
    ct.save_dir_name = 'xe_history_hubble' + '_G_' + str_gamma

    R = ct.run_cosmotherm_parallel(**args)
print(R)
