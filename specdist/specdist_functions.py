from .utils import *
from .cosmology import *


def redshift_z_mu(cosmo):
    #see eq. 4.47 of https://physique.cuso.ch/fileadmin/physique/document/2014_Chluba_notes.pdf
    #this assumes N_eff = 3.046
    #this is only the double compton thermalization redshift
    #return 1.98e6*(cosmo.omega_b/0.022)**-(2./5.)*((1.-cosmo.Yp/2.)/0.88)**-(2./5.)*(cosmo.T_cmb/2.725)**(1./5.)
    return 1.98e6

def visibility_J_bb(z,cosmo):
    #eq. 4.46 of https://physique.cuso.ch/fileadmin/physique/document/2014_Chluba_notes.pdf
    #this is assuming DC only
    #z = np.asarray(z)
    try:
        result = np.exp(-(z/redshift_z_mu(cosmo))**(5./2.))
    except:
        result = 0.
    if math.isnan(result):
        result = 0.
    return result

def visibility_J_bb_star(z,cosmo):
    #see eq. 13 of https://arxiv.org/pdf/1506.06582.pdf
    try:
        result = 0.983*visibility_J_bb(z,cosmo)*(1.-0.0381*(z/redshift_z_mu(cosmo))**2.29)
    except:
        result = 0.
    if math.isnan(result):
        result = 0.
    return result

def visibility_J_y(z,cosmo):
    #see eq. 5 of https://arxiv.org/pdf/1304.6120.pdf
    #z = np.asarray(z)
    result = (1.+((1.+z)/6e4)**2.58)**-1.
    if math.isnan(result):
        result = 0.
    return result

def visibility_J_mu(z,cosmo):
    #see eq. 5 of https://arxiv.org/pdf/1304.6120.pdf
    try:
        result = 1.-np.exp(-((1.+z)/5.8e4)**1.88)
    except:
        result = 0.
    if math.isnan(result):
        result = 0.
    return result

def visibility_J_T(z,cosmo):
    #see eq. 5 of https://arxiv.org/pdf/1304.6120.pdf
    result = 1.-visibility_J_bb_star(z,cosmo)
    if math.isnan(result):
        result = 0.
    return result

def critical_frequency_x_c_br(z):
    #eq. 4.39 of https://physique.cuso.ch/fileadmin/physique/document/2014_Chluba_notes.pdf
    #assumes Itoh et al BR treatment
    return 1.23e-3*((1.+z)/2e6)**-0.672

def critical_frequency_x_c_dc(z):
    #eq. 4.38 of https://physique.cuso.ch/fileadmin/physique/document/2014_Chluba_notes.pdf
    #assumes DC Gaunt factors are negligible
    return 8.60e-3*((1.+z)/2e6)**0.5

def critical_frequency_x_c(z):
    return np.sqrt(critical_frequency_x_c_br(z)**2.+critical_frequency_x_c_dc(z)**2.)


def mu_from_energy_release_history(energy_release_history_dlnrho_dt,cosmo,**kwargs):
    def integrand(ln1pz,*args):
        z = np.exp(ln1pz)-1.
        J_bb = visibility_J_bb(z,args[0])
        J_mu = visibility_J_mu(z,args[0])
        dt_dln1pz = -1./cosmo.E(z)/args[0].H0()
        dlnrho_dln1pz = energy_release_history_dlnrho_dt(z,args[0],**args[1])*dt_dln1pz
        result = 3./kappa_c*J_bb*J_mu*dlnrho_dln1pz
        return result
    #trapezoidal rule
    nz = int(50)
    ln1pz_array = np.linspace((np.log(1.+cosmo.z_start)),(np.log(1.+cosmo.z_end)),nz)
    Ip = []
    int_array_xp = []
    a_args = (cosmo,kwargs)
    for p in ln1pz_array:
        int_p = integrand(p,*a_args)
        int_array_xp.append(int_p)
    int_array_xp=np.asarray(int_array_xp)
    Ip = np.trapz(int_array_xp,ln1pz_array)
    result = (Ip,0.)
    ####end trapezoidal rule
    #result =  quad(integrand,np.log(1.+cosmo.z_start),np.log(1.+cosmo.z_end), args=(cosmo,kwargs))
    r_dict = {}
    r_dict['value']=result[0]
    r_dict['err'] = result[1]
    return r_dict



def y_from_energy_release_history(energy_release_history_dlnrho_dt,cosmo,**kwargs):
    def integrand(ln1pz,*args):
        z = np.exp(ln1pz)-1.
        J_bb = visibility_J_bb(z,args[0])
        J_y = visibility_J_y(z,args[0])
        dt_dln1pz = -1./cosmo.E(z)/args[0].H0()
        dlnrho_dln1pz = energy_release_history_dlnrho_dt(z,args[0],**args[1])*dt_dln1pz
        result = J_bb*J_y*dlnrho_dln1pz/4.
        return result
    #trapezoidal rule
    nz = int(50)
    ln1pz_array = np.linspace((np.log(1.+cosmo.z_start)),(np.log(1.+cosmo.z_end)),nz)
    Ip = []
    int_array_xp = []
    a_args = (cosmo,kwargs)
    for p in ln1pz_array:
        int_p = integrand(p,*a_args)
        int_array_xp.append(int_p)
    int_array_xp=np.asarray(int_array_xp)
    Ip = np.trapz(int_array_xp,ln1pz_array)
    result = (Ip,0.)
    ####end trapezoidal rule
    #result =  quad(integrand,np.log(1.+cosmo.z_start),np.log(1.+cosmo.z_end), args=(cosmo,kwargs))
    r_dict = {}
    r_dict['value']=result[0]
    r_dict['err'] = result[1]
    return r_dict

def Drho_rho_y_from_energy_release_history(energy_release_history_dlnrho_dt,cosmo,**kwargs):
    return y_from_energy_release_history(energy_release_history_dlnrho_dt,cosmo,**kwargs)['value']*4.

def Drho_rho_mu_from_energy_release_history(energy_release_history_dlnrho_dt,cosmo,**kwargs):
    return mu_from_energy_release_history(energy_release_history_dlnrho_dt,cosmo,**kwargs)['value']/(3./kappa_c)

def Drho_rho_tot_from_energy_release_history(energy_release_history_dlnrho_dt,cosmo,**kwargs):
    return Drho_rho_y_from_energy_release_history(energy_release_history_dlnrho_dt,cosmo,**kwargs)+Drho_rho_mu_from_energy_release_history(energy_release_history_dlnrho_dt,cosmo,**kwargs)


def DN_N_from_entropy_production_history(entropy_production_history_dlnN_dt,cosmo,**kwargs):
    def integrand(ln1pz,*args):
        z = np.exp(ln1pz)-1.
        # J_bb = visibility_J_bb(z,args[0])
        # J_y = visibility_J_y(z,args[0])
        dt_dln1pz = -1./cosmo.E(z)/args[0].H0()
        dlnN_dln1pz = entropy_production_history_dlnN_dt(z,args[0],**args[1])*dt_dln1pz
        result = dlnN_dln1pz
        return result
    #trapezoidal rule
    nz = int(50)
    ln1pz_array = np.linspace((np.log(1.+cosmo.z_start)),(np.log(1.+cosmo.z_end)),nz)
    Ip = []
    int_array_xp = []
    a_args = (cosmo,kwargs)
    for p in ln1pz_array:
        int_p = integrand(p,*a_args)
        int_array_xp.append(int_p)
    int_array_xp=np.asarray(int_array_xp)
    Ip = np.trapz(int_array_xp,ln1pz_array)
    result = (Ip,0.)
    ####end trapezoidal rule
    #result =  quad(integrand,np.log(1.+cosmo.z_start),np.log(1.+cosmo.z_end), args=(cosmo,kwargs))
    r_dict = {}
    r_dict['value']=result[0]
    r_dict['err'] = result[1]
    return r_dict
