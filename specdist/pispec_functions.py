from .utils import *
from .cosmology import *
from .cosmotherm_wrapper import *
from .specdist_functions import *

class dm_particle:
    f_gamma = 2
    f_dm = 1.
    Gamma_inj = 1e-17
    x_0 = 1e-3

def x_inj(dm_particle,z):
    return dm_particle.x_0/(1.+z)


def f_inj(dm_particle,cosmo):
    return 4.849e3*dm_particle.f_dm*dm_particle.f_gamma/dm_particle.x_0*(cosmo.omega_cdm/0.12)*(cosmo.T_cmb/2.726)**-4



def injection_redshift_zX(cosmo,cosmotherm):
    # set_cosmo_to_CT_cosmo_params(cosmo,cosmotherm)

    def f(lnz):
        z = np.exp(lnz)
        return cosmotherm.ct_Gamma_dec*cosmo.t_H0_in_s()-cosmo.E(z)-cosmo.dE_dz(z)*(1.+z)
    root = np.exp(optimize.brentq(f, np.log(cosmotherm.ct_zstart), np.log(cosmotherm.ct_zend)))
    return root

def mu_instantaneous_injection(zi,cosmo,cosmotherm,dm_particle):
    x_0 = 4./3./a_rho
    x_i = x_inj(dm_particle,0.)

    DN_N = 1.
    x_c = critical_frequency_x_c(zi)
    P_s = np.exp(-x_c/x_i)

    term_1 = 0.#3./kappa_c*cosmotherm.ct_Drho_rho_dec
    #term_2 = 3.*a_rho/kappa_c*(x_i-x_0)*visibility_J_bb_star(zi,cosmo)*DN_N
    term_2 = 3.*a_rho/kappa_c*(x_i-x_0*P_s)*visibility_J_bb_star(zi,cosmo)*DN_N

    return term_1 + term_2

def dI_dln1pz_of_z(z,cosmo,cosmotherm):
    # set_cosmo_to_CT_cosmo_params(cosmo,cosmotherm)
    delta_t = cosmo.t_of_z_in_s(z)['value'] - cosmo.t_of_z_in_s(cosmo.z_start)['value']
    return np.exp(-cosmotherm.ct_Gamma_dec*delta_t)/(1.+z)/cosmo.E(z)


def I(cosmo,cosmotherm):
    def integrand(ln1pz,cosmo,cosmotherm):
        z = np.exp(ln1pz)-1.
        return  dI_dln1pz_of_z(z,cosmo,cosmotherm)
    #zend = np.maximum(cosmo.z_end,cosmo.z_of_t(10.*1./cosmotherm.ct_Gamma_dec))
    zend = cosmo.z_end
    # print(zend)
    result =  quad(integrand,np.log(1.+zend),np.log(1.+cosmo.z_start), args=(cosmo,cosmotherm))
    r_dict = {}
    r_dict['value']=result[0]
    r_dict['err'] = result[1]
    return r_dict

def Drho_rho_inj_at_z_normalized(z,cosmo,cosmotherm):
    # print(cosmotherm.ct_Gamma_dec)
    # print(cosmo.Omega_m())
    # print(cosmo.Omega_rel()*cosmo.h**2.)
    z = np.asarray(z)
    if len(z) ==1:
        return dI_dln1pz_of_z(z,cosmo,cosmotherm)/I(cosmo,cosmotherm)['value']
    else:
        denominator = I(cosmo,cosmotherm)['value']
        numerator = []
        for zp in z:
            numerator.append(dI_dln1pz_of_z(zp,cosmo,cosmotherm))
        numerator = np.asarray(numerator)
        return numerator/denominator

def pi_energy_release_history_dlnrho_dt(z,cosmo,**kwargs):
    ct = kwargs['cosmotherm']
    X_dm = kwargs['dm_particle']
    delta_t = cosmo.t_of_z_in_s(z)['value'] - cosmo.t_of_z_in_s(cosmo.z_start)['value']
    #note: this is independent of x_inj
    return G2/G3*f_inj(X_dm,cosmo)*X_dm.Gamma_inj*x_inj(X_dm,z)*np.exp(-X_dm.Gamma_inj*delta_t)


def set_dm_params_to_CT_pi_params(dm_particle,cosmotherm):
    dm_particle.Gamma_inj = cosmotherm.ct_Gamma_dec
    dm_particle.x_0 = cosmotherm.ct_x_dec

def get_fdm_from_Drho_rho_tot(Drho_rho_tot,cosmo,cosmotherm,dm_particle):
    dm_particle.f_dm = 1.
    dict = {}
    dict['cosmotherm']=cosmotherm
    dict['dm_particle']=dm_particle
    return Drho_rho_tot/Drho_rho_tot_from_energy_release_history(pi_energy_release_history_dlnrho_dt,cosmo,**dict)
