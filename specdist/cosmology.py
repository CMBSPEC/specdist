from .utils import *

class cosmo:
    T_cmb = 2.726 #Kelvins
    omega_cdm = 0.120 #Planck 2018
    omega_b = 0.0224 #Planck 2018
    h = 0.71
    z_start = 5e6
    z_end = 1e-2
    N_eff = 3.046
    Yp = 0.24

    def H0(self):
        return 100.*self.h*km_over_Mpc

    def Omega_cdm(self):
        return self.omega_cdm/self.h**2
    def Omega_b(self):
        return self.omega_b/self.h**2
    def Omega_m(self):
        return self.Omega_cdm()+self.Omega_b()
    def Omega_g(self):
        return omega_gamma(self.T_cmb)/self.h**2
    def Omega_rel(self):
        a=np.pi**2/15.0*kb*self.T_cmb*np.power(kb*self.T_cmb/(hplanck*clight/2./np.pi), 3)/np.power(clight, 2)
        b=3.0*np.power(self.H0(), 2)/(8.0*np.pi*G_newton)
        # print(a,b,kb,hplanck,clight,self.H0())
        return  a/b*(1.0+self.N_eff*(7.0/8.0)*np.power(4.0/11.0, 4.0/3.0))

    def z_eq(self):
        return self.Omega_m()/self.Omega_rel()-1.

    def Omega_L(self):
        return 1.-self.Omega_m() -self.Omega_rel()
    def E(self,z):
        return np.sqrt(self.Omega_m()*(1.+z)**3.+self.Omega_rel()*(1.+z)**4.+self.Omega_L())
    def dE_dz(self,z):
        denominator = 2.*self.E(z)
        numerator = 3.*self.Omega_m()*(1.+z)**2.+4.*self.Omega_rel()*(1.+z)**3.
        return numerator/denominator

    def t_H0_in_s(self):
        return 1./self.H0()

    def dt_dz_of_z_in_s(self,z):
        return -self.t_H0_in_s()*1./self.E(z)/(1.+z)

    def t_of_z_in_s(self,z):
        def integrand(ln1pz,self):
            zp = np.exp(ln1pz)-1.
            dzdln1pz= 1.+zp
            return self.dt_dz_of_z_in_s(zp)*dzdln1pz
        result =  quad(integrand,np.log(1.+1e14),np.log(1.+z), args=self)
        r_dict = {}
        r_dict['value']=result[0]
        r_dict['err'] = result[1]
        return r_dict

    def z_of_t(self,t):
        if t<=self.t_of_z_in_s(cosmo.z_start)['value']:
            return cosmo.z_start
        if t>=self.t_of_z_in_s(cosmo.z_end)['value']:
            return cosmo.z_end
        def f(ln1pz):
            z = np.exp(ln1pz)-1.
            return self.t_of_z_in_s(z)['value']-t
        root = np.exp(optimize.brentq(f, np.log(1.+self.z_start), np.log(1.+self.z_end)))-1.
        return root




def set_cosmo_to_CT_cosmo_params(cosmo,cosmotherm):
    cosmo.T_cmb = cosmotherm.ct_T0
    cosmo.h = cosmotherm.ct_h
    cosmo.omega_b = cosmotherm.ct_Omega_b*cosmo.h**2.
    cosmo.omega_cdm = (cosmotherm.ct_Omega_m - cosmotherm.ct_Omega_b)*cosmo.h**2.
    cosmo.z_start = cosmotherm.ct_zstart
    cosmo.z_end = cosmotherm.ct_zend
    cosmo.N_eff = cosmotherm.ct_N_eff
    cosmo.Yp = cosmotherm.ct_Yp



def rho_gamma_in_GeV_per_cm3(T):
    return 8.*np.pi**5.*(kb*T)**4/15./clight**3./hplanck**3./GeV_over_kg*1e-6/clight**2


def omega_gamma(T):
    return rho_gamma_in_GeV_per_cm3(T)/rho_crit_over_h2_in_GeV_per_cm3
