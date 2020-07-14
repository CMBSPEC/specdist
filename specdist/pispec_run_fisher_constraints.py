from .utils import *
from .cosmology import *
from .cosmotherm_wrapper import *
from .specdist_functions import *
from .pispec_functions import *
from .interpolate_ct_spectra import *
from .standard_mu_and_y_distortions import *
from .firas_data import *
from scipy.linalg import cholesky, LinAlgError
from scipy.linalg import block_diag

def pi_run_fisher_constraints(Gamma_values,xi_array,sd_lib,*args,**kwargs):
    add_edges = kwargs.get('add_edges', 'no')
    #
    # photon_injection_case = 'lyc'
    # sd_lib = sd.specdist_ct_spectra_lib()
    # sd.load_ct_spectra_lib(photon_injection_case,sd_lib)
    # sd_lib_lyc = sd_lib
    #
    # photon_injection_case = 'bare'
    # sd_lib = sd.specdist_ct_spectra_lib()
    # sd.load_ct_spectra_lib(photon_injection_case,sd_lib)
    # sd_lib_bare = sd_lib
    #
    # photon_injection_case = 'lyc_reio'
    # sd_lib = sd.specdist_ct_spectra_lib()
    # sd.load_ct_spectra_lib(photon_injection_case,sd_lib)
    # sd_lib_lyc_reio = sd_lib


    firas = kwargs['firas']
    edges = kwargs['edges']


    T_pivot = 2.725
    GetFirasXTatTpivot(firas,T_pivot)

    # photon_injection_case = 'lyc_reio'
    # add_edges = 'no'

    # sd_lib = sd_lib_lyc_reio
    #theory
    #parameters
    id_param = 0
    theta_parameters = []
    param_delta_T = 1.
    id_param_delta_T = id_param
    id_param += 1
    theta_parameters.append(param_delta_T)
    param_G0 = 1.
    id_param_G0 = id_param
    id_param += 1
    theta_parameters.append(param_G0)
    param_finj_over_fct = 1.
    id_param_finj_over_fct = id_param
    id_param += 1
    theta_parameters.append(param_finj_over_fct)

    M_parameters = len(theta_parameters)



    firas_delta_monopole = []
    firas_galaxy = []
    theory_photon_injection_spectrum = []

    for i in range(len(firas.firas_x)):
        nu = firas.firas_nu[i]
        firas_delta_monopole.append(dB_nu_dT_at_T(nu,T_pivot))
        firas_galaxy.append(firas.firas_Gnu[i])

    covmat_data = firas.firas_covmat

    #Add edges_like point
    if add_edges=='yes':
        #edges = sd.edges()
        I_edges = B_nu_of_T(edges.edges_nu,T_pivot)
        firas.firas_x = np.insert(firas.firas_x,0,edges.edges_x)
        firas_galaxy = np.insert(firas_galaxy,0,0.)
        firas_delta_monopole.append(dB_nu_dT_at_T(edges.edges_nu,T_pivot))
        cov_edges = [[(2.*I_edges)**2.]]
        covmat_data = block_diag(cov_edges, covmat_data)

    #covmat_data = np.diag(covmat_data)
    inv_covmat_data = np.linalg.inv(covmat_data)
    det_covmat_data = np.linalg.det(covmat_data)
    try:
        cholesky(covmat_data)
    except LinAlgError:
        raise ValueError("Covariance is not SPD!")


    x_asked = firas.firas_x
    # gamma_asked = 1e-8
    # Nx = 400
    # xi_array = np.logspace(-6,6,Nx)

    #Gamma_values = Gamma_inj_asked
    #gamma_labels = [r'$10^{-8}$',r'$10^{-9}$',r'$10^{-10}$',r'$10^{-11}$',r'$10^{-12}$',r'$10^{-13}$',r'$10^{-14}$',r'$10^{-15}$',r'$10^{-16}$',r'$10^{-17}$']


    f_dm_fisher = {}
    f_dm_fisher['curves'] = []
    f_dm_fisher['Gamma_inj'] = []

    for gamma_asked in Gamma_values:
        gamma = gamma_asked
        mu_high_z = {}
        mu_high_z['x'] = xi_array
        mui = []
        finji = []
        for xinj_asked in xi_array:
            S = GetSpectra(gamma_asked,xinj_asked,x_asked,sd_lib)
            theory_photon_injection_spectrum = S['DI']*1e-6
            #print(S['DI'])


            mu_parameters = []
            mu_parameters.append(firas_delta_monopole)
            mu_parameters.append(firas_galaxy)
            mu_parameters.append(theory_photon_injection_spectrum)
            mu_parameters = np.asarray(mu_parameters)


            dmu_parameters = []
            for m in range(M_parameters):
                dmu_parameters.append(np.asarray(mu_parameters[m]))


            fisher_F = np.zeros((M_parameters,M_parameters))
            for a in range(M_parameters):
                for b in range(M_parameters):
                    fisher_M_ab = np.outer(dmu_parameters[a],dmu_parameters[b]) + np.outer(dmu_parameters[b],dmu_parameters[a])
                    fisher_F[a][b]=0.5*np.trace(np.matmul(inv_covmat_data,fisher_M_ab))

            inverse_fisher_F = np.linalg.inv(fisher_F)
            fisher_sigmas = []
            for m in range(M_parameters):
                fisher_sigmas.append(np.sqrt(inverse_fisher_F[m][m]))
            fisher_sigmas = np.asarray(fisher_sigmas)

            #mui.append(fisher_sigmas[id_param_finj_over_fct])
            r_finj_over_fct =  2.*fisher_sigmas[id_param_finj_over_fct]
            f_ct =  S['finj']
            f_dm = 1./1.3098e4*r_finj_over_fct*f_ct*xinj_asked
            mui.append(f_dm)
            finji.append(r_finj_over_fct*f_ct)
        mu_high_z['fdm'] = np.asarray(mui)
        mu_high_z['finj'] = np.asarray(finji)
        f_dm_fisher['curves'].append(mu_high_z)
        f_dm_fisher['Gamma_inj'].append(gamma)


    return f_dm_fisher
    #print(f_dm_fisher['curves'][0]['fdm'])
