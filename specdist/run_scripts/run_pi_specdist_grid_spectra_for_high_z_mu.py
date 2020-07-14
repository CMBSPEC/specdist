import specdist as pi
import numpy as np

# load photon injection and cosnotherm modules
ct = pi.cosmotherm()
cosmo = pi.cosmo()
X_dm = pi.dm_particle()

# set relevant parameter values
ct.ct_Drho_rho_dec = 1e-100
ct.ct_h = 0.70
ct.ct_Omega_b = 0.0457
ct.ct_Omega_m = 0.30
ct.ct_emission_absorption_mode = 4
ct.ct_npts = 5000
ct.ct_zend = 1e-2
pi.set_dm_params_to_CT_pi_params(X_dm,ct)
pi.set_cosmo_to_CT_cosmo_params(cosmo,ct)


xi_array = np.logspace(np.log10(1e1),np.log10(1.e8),50)

#Gamma_X = 1e-6 #np.logspace(-9,-6,4)

for Gamma_X in np.logspace(-9,-6,4):
    str_gamma = str("%.3e"%Gamma_X)

    ct.ct_Gamma_dec = Gamma_X



    args = {}
    p_name = 'photon injection x_dec'
    p_array = xi_array
    args['param_values_array'] = p_array
    args['param_name'] = p_name
    args['save_spectra'] = 'yes'
    ct.save_dir_name = 'spectra_for_mu_fit_hubble' + '_G_' + str_gamma

    R = ct.run_cosmotherm_parallel(**args)
