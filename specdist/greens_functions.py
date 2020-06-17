from .utils import *
from .specdist_functions import *
from .standard_mu_and_y_distortions import *

def greens_functions_G_in_MJy_per_sr(x,z,cosmo):
    #eq 6 of https://arxiv.org/pdf/1304.6120.pdf
    mu_part =  3./kappa_c*visibility_J_bb_star(z,cosmo)*visibility_J_mu(z,cosmo)*M(x)
    y_part = 1./4.*visibility_J_y(z,cosmo)*Y_sz(x)
    T_part = (1.-visibility_J_bb_star(z,cosmo))/4.*G_bb(x)
    return DI_normalization_in_MJy_per_sr(x,cosmo)*(mu_part+y_part+T_part)


def greens_functions_DI_from_energy_release_history_in_MJy_per_sr(X,energy_release_history_dlnrho_dt,cosmo,**kwargs):
    def integrand(ln1pz,*args):
        z = np.exp(ln1pz)-1.
        dt_dln1pz = -1./cosmo.E(z)/args[0].H0()
        dlnrho_dln1pz = energy_release_history_dlnrho_dt(z,args[0],**args[1])*dt_dln1pz
        x = args[2]
        result = greens_functions_G_in_MJy_per_sr(x,z,cosmo)*dlnrho_dln1pz
        return result
    x = np.asarray(X)
    dist = []
    dist_err = []
    #trapezoidal rule
    nz = int(50)
    #ln1pz_array = np.logspace(np.log10(np.log(1.+cosmo.z_start)),np.log10(np.log(1.+cosmo.z_end)),nz)
    ln1pz_array = np.linspace((np.log(1.+cosmo.z_start)),(np.log(1.+cosmo.z_end)),nz)
    Ip = []

    for xp in x:
        a_args = (cosmo,kwargs,xp)
        int_array_xp = []
        for p in ln1pz_array:
            int_p = integrand(p,*a_args)
            int_array_xp.append(int_p)
        int_array_xp=np.asarray(int_array_xp)
        Ip.append(np.trapz(int_array_xp,ln1pz_array))
    Ip = np.asarray(Ip)
    return (Ip,np.zeros(len(Ip)))
    try:
        for xp in x:
            result =  quad(integrand,np.log(1.+cosmo.z_start),np.log(1.+cosmo.z_end), args=(cosmo,kwargs,xp))
            dist.append(result[0])
            dist_err.append(result[1])
    except:
        result =  quad(integrand,np.log(1.+cosmo.z_start),np.log(1.+cosmo.z_end), args=(cosmo,kwargs,x))
        dist.append(result[0])
        dist_err.append(result[1])

    return (np.asarray(dist),np.asarray(dist_err))
