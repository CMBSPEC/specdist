from .utils import *

class cosmo:
    T_cmb = 2.726 #Kelvins
    omega_cdm = 0.120 #Planck 2018
    omega_b = 0.0224 #Planck 2018
    h = 0.71

    def Omega_cdm(self):
        return self.omega_cdm/self.h**2
    def Omega_b(self):
        return self.omega_b/self.h**2
    def Omega_m(self):
        return self.Omega_cdm()+self.Omega_b()
    def Omega_g(self):
        return omega_gamma(self.T_cmb)/self.h**2
    def Omega_L(self):
        return 1.-self.Omega_m()-self.Omega_g()
    def E(self,z):
        return np.sqrt(self.Omega_m()*(1.+z)**3.+self.Omega_g()*(1.+z)**4.+self.Omega_L())
    def dE_dz(self,z):
        denominator = 2.*self.E(z)
        numerator = 3.*self.Omega_m()*(1.+z)**2.+4.*self.Omega_g()*(1.+z)**3.
        return numerator/denominator

    def t_H0_in_s(self):
        return 1./km_over_Mpc/self.h*1e-2


def set_cosmo_to_CT_cosmo_params(cosmo,cosmotherm):
    cosmo.T_cmb = cosmotherm.ct_T0
    cosmo.h = cosmotherm.ct_h
    cosmo.omega_b = cosmotherm.ct_Omega_b*cosmo.h**2.
    cosmo.omega_cdm = (cosmotherm.ct_Omega_m - cosmotherm.ct_Omega_b)*cosmo.h**2.





rho_crit_over_h2_in_GeV_per_cm3 = 1.0537e-5



def rho_gamma_in_GeV_per_cm3(T):
    return 8.*np.pi**5.*(kb*T)**4/15./clight**3./hplanck**3./GeV_over_kg*1e-6/clight**2


def omega_gamma(T):
    return rho_gamma_in_GeV_per_cm3(T)/rho_crit_over_h2_in_GeV_per_cm3
