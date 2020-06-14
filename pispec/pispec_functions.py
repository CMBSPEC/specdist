from .utils import *
from .cosmology import *

def injection_redshift_zX(cosmo,cosmotherm):
    set_cosmo_to_CT_cosmo_params(cosmo,cosmotherm)

    def f(lnz):
        z = np.exp(lnz)
        return cosmotherm.ct_Gamma_dec*cosmo.t_H0_in_s()-cosmo.E(z)-cosmo.dE_dz(z)*(1.+z)
    root = np.exp(optimize.brentq(f, np.log(cosmotherm.ct_zstart), np.log(cosmotherm.ct_zend)))
    return root

def mu_instantaneous_injection(zi,cosmo,cosmotherm):
    alpha_rho = 4.*0.0925
    x_0 = 4./3./alpha_rho
    x_i = cosmotherm.ct_x_dec/(1.+zi)

    DN_N = 2
    x_c = np.sqrt(1.)
    P_s = np.exp(-x_c/x_i)

    term_1 = 3./kappa_c*cosmotherm.ct_Drho_rho_dec
    term_2 = 3.*alpha_rho/kappa_c*(x_i-x_0*P_s)*DN_N

    return term_1 + term_2
