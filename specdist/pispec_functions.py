from .utils import *
from .cosmology import *
from .cosmotherm_wrapper import *
from .specdist_functions import *


#Particle mass in eV, m/1eV = 4.698*10^-4 x_inj
class dm_particle:
    f_gamma = 2
    f_dm = 1.
    Gamma_inj = 1e-17
    x_0 = 1e-3

def x_inj(dm_particle,z):
    return dm_particle.x_0/(1.+z)


def f_inj(dm_particle,cosmo):
    #factor 2 when f_gamma = 2
    return 1./2.*1.3098e4*dm_particle.f_dm*dm_particle.f_gamma/dm_particle.x_0*(cosmo.omega_cdm/0.12)*(cosmo.T_cmb/2.726)**-4



def get_f_inj_from_Drho_rho(Drho_rho_tot,cosmo,cosmotherm,dm_particle):
    dm_particle.f_dm = get_fdm_from_Drho_rho(Drho_rho_tot,cosmo,cosmotherm,dm_particle)['tot']
    #factor 2 when f_gamma = 2
    return 1./2.*1.3098e4*dm_particle.f_dm*dm_particle.f_gamma/dm_particle.x_0*(cosmo.omega_cdm/0.12)*(cosmo.T_cmb/2.726)**-4



def injection_redshift_zX(gamma_inj,cosmo,cosmotherm):
    # set_cosmo_to_CT_cosmo_params(cosmo,cosmotherm)

    def f(lnz):
        z = np.exp(lnz)
        return gamma_inj*cosmo.t_H0_in_s()-cosmo.E(z)-cosmo.dE_dz(z)*(1.+z)
    root = np.exp(optimize.brentq(f, np.log(cosmotherm.ct_zstart), np.log(cosmotherm.ct_zend)))
    return root

def find_Gamma_inj_for_injection_redshift_zX(zX,cosmo,cosmotherm):
    def f(lng):
        g = np.exp(lng)
        return g*cosmo.t_H0_in_s()-cosmo.E(zX)-cosmo.dE_dz(zX)*(1.+zX)
    root = np.exp(optimize.brentq(f, np.log(1e-20), np.log(1e-2)))
    return root


def high_redshift_f_dm_limit(mu_lim,cosmo,cosmotherm,dm_particle,*args,**kwargs):
    N_int = kwargs.get('N_int', 50)
    # def f(ln_fdm):
    #     f_dm = np.exp(ln_fdm)
    #     dm_particle.f_dm = f_dm
    #     mui = np.abs(mu_continuous_injection(cosmo,cosmotherm,dm_particle,N_int = N_int)['value'])
    #     return mui-mu_lim
    # root = np.exp(optimize.brentq(f, np.log(1e-100), np.log(1e100)))
    # return mu_lim

    dm_particle.f_dm = 1.
    mui = mu_continuous_injection(cosmo,cosmotherm,dm_particle,N_int = N_int)['value']
    #return mui-mu_lim
    #root = np.exp(optimize.brentq(f, np.log(1e-100), np.log(1e100)))
    return mu_lim/mui



def mu_instantaneous_injection(zi,cosmo,cosmotherm,dm_particle):
    x_0 = 4./3./a_rho
    x_i = x_inj(dm_particle,zi)

    DN_N = 1.
    x_c = critical_frequency_x_c(zi)
    P_s = np.exp(-x_c/x_i)

    term_1 = 0.#3./kappa_c*cosmotherm.ct_Drho_rho_dec
    #term_2 = 3.*a_rho/kappa_c*(x_i-x_0)*visibility_J_bb_star(zi,cosmo)*DN_N
    term_2 = 3.*a_rho/kappa_c*(x_i-x_0*P_s)*visibility_J_bb_star(zi,cosmo)*DN_N
    return term_1 + term_2




def dmu_dt_continuous_injection(zi,cosmo,**kwargs):
    ct = kwargs['cosmotherm']
    X_dm = kwargs['dm_particle']
    x_0 = 4./3./a_rho
    x_i = x_inj(X_dm,zi)
    DN_N = pi_entropy_production_history_dlnN_dt(zi,cosmo,**kwargs)
    x_c = critical_frequency_x_c(zi)
    P_s = np.exp(-x_c/x_i)

    term_1 = 0.#3./kappa_c*cosmotherm.ct_Drho_rho_dec
    #term_2 = 3.*a_rho/kappa_c*(x_i-x_0)*visibility_J_bb_star(zi,cosmo)*DN_N
    term_2 = 3.*a_rho/kappa_c*(x_i-x_0*P_s)*visibility_J_bb_star(zi,cosmo)*DN_N
    return term_1 + term_2

def mu_continuous_injection(cosmo,cosmotherm,dm_particle,*args,**kwargs):
    N_int = kwargs.get('N_int', 50)
    dict = {}
    dict['cosmotherm']=cosmotherm
    dict['dm_particle']=dm_particle
    def integrand(ln1pz,*args):
        z = np.exp(ln1pz)-1.
        dt_dln1pz = -1./cosmo.E(z)/args[0].H0()
        dmu_dln1pz = dmu_dt_continuous_injection(z,args[0],**args[1])*dt_dln1pz
        result = dmu_dln1pz
        return result

    #trapezoidal rule
    nz = int(N_int)
    zend_gamma = cosmo.z_of_t(10.*1./dm_particle.Gamma_inj)
    zend = max(cosmo.z_end,zend_gamma)
    ln1pz_array = np.linspace((np.log(1.+cosmo.z_start)),(np.log(1.+zend)),nz)
    Ip = []
    int_array_xp = []
    a_args = (cosmo,dict)
    for p in ln1pz_array:
        int_p = integrand(p,*a_args)
        int_array_xp.append(int_p)
    int_array_xp=np.asarray(int_array_xp)
    Ip = np.trapz(int_array_xp,ln1pz_array)
    result = (Ip,0.)
    ####end trapezoidal rule
    #result =  quad(integrand,np.log(1.+cosmo.z_start),np.log(1.+cosmo.z_end), args=(cosmo,dict))
    r_dict = {}
    r_dict['value']=result[0]
    r_dict['err'] = result[1]
    return r_dict






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
    tz = cosmo.t_of_z_in_s(z)['value']
    delta_t = tz - cosmo.t_of_z_in_s(cosmo.z_start)['value']
    if X_dm.Gamma_inj*delta_t>100.:
        return 0.
    else:
        #note: this is independent of x_inj
        return G2/G3*f_inj(X_dm,cosmo)*X_dm.Gamma_inj*x_inj(X_dm,z)*np.exp(-X_dm.Gamma_inj*delta_t)


def set_dm_params_to_CT_pi_params(dm_particle,cosmotherm):
    dm_particle.Gamma_inj = cosmotherm.ct_Gamma_dec
    dm_particle.x_0 = cosmotherm.ct_x_dec

def get_fdm_from_Drho_rho(Drho_rho_tot,cosmo,cosmotherm,dm_particle):
    dm_particle.f_dm = 1.
    dict = {}
    dict['cosmotherm']=cosmotherm
    dict['dm_particle']=dm_particle
    numerator = Drho_rho_tot
    denominator_y = Drho_rho_y_from_energy_release_history(pi_energy_release_history_dlnrho_dt,cosmo,**dict)
    denominator_mu = Drho_rho_mu_from_energy_release_history(pi_energy_release_history_dlnrho_dt,cosmo,**dict)
    r_dict = {}
    try:
        r_dict['y'] = numerator/denominator_y
    except:
        r_dict['y'] = 0.
    try:
        r_dict['mu'] = numerator/denominator_mu
    except:
        r_dict['mu'] = 0.
    try:
        r_dict['tot'] = numerator/(denominator_mu+denominator_y)
    except:
        r_dict['tot'] = 0.
    return r_dict

def pi_entropy_production_history_dlnN_dt(z,cosmo,**kwargs):
    ct = kwargs['cosmotherm']
    X_dm = kwargs['dm_particle']
    delta_t = cosmo.t_of_z_in_s(z)['value'] - cosmo.t_of_z_in_s(cosmo.z_start)['value']
    #note: this is *dependent* of x_inj
    return f_inj(X_dm,cosmo)*X_dm.Gamma_inj*np.exp(-X_dm.Gamma_inj*delta_t)

def get_fdm_from_mu_continuous_injection(mu_limit,cosmo,cosmotherm,dm_particle):
    dm_particle.f_dm = 1.
    numerator = mu_limit
    denominator = mu_continuous_injection(cosmo,cosmotherm,dm_particle)['value']
    try:
        result = numerator/denominator
    except:
        result = 0.
    return result
