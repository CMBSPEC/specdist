from .utils import *





#Black body spectrum (in MJy/sr)
def B_nu_of_T(NU,T):
    return (2.*hplanck/clight**2.)*NU**3./(np.exp(hplanck*NU/kb/T)-1.)*1.e20

#Derivative of black body spectrum (in MJy/sr/K)
def dB_nu_dT_at_T(NU,T):
    return B_nu_of_T(NU,T)/T*(hplanck*NU/kb/T)*np.exp(hplanck*NU/kb/T)/(np.exp(hplanck*NU/kb/T)-1.)

#MU distortion
def dS_dMU(NU,T):
    x = (hplanck*NU/kb/T)
    return -T/x*dB_nu_dT_at_T(NU,T)

#Y distortion
def dS_dY(NU,T):
    x = (hplanck*NU/kb/T)
    return T*(x/np.tanh(x/2.)-4.)*dB_nu_dT_at_T(NU,T)


def GetMuSpecDistAtTandX(mu,T,X):
    x = np.asarray(X)
    dist = []
    try:
        for xp in x:
            nu = kb*T/hplanck*xp
            dist.append(mu*dS_dMU(nu,T))
    except:
        nu = kb*T/hplanck*x
        dist.append(mu*dS_dMU(nu,T))
    return np.asarray(dist)

def GetMuSpecDistAtTandX_chluba(mu,T,X):
    x = np.asarray(X)
    dist = []
    try:
        for xp in x:
            nu = kb*T/hplanck*xp
            dist.append(mu*(-dS_dMU(nu,T)*xp*(1./beta_mu-1./xp)))
    except:
        nu = kb*T/hplanck*x
        dist.append(mu*(-dS_dMU(nu,T)*x*(1./beta_mu-1./x)))
    return np.asarray(dist)

def GetYSpecDistAtTandX(y,T,X):
    x = np.asarray(X)
    dist = []
    try:
        for xp in x:
            nu = kb*T/hplanck*xp
            dist.append(y*dS_dY(nu,T))
    except:
        nu = kb*T/hplanck*x
        dist.append(y*dS_dY(nu,T))
    return np.asarray(dist)

def G_bb(x):
    result =  x*np.exp(x)/(np.exp(x)-1.)**2.
    if math.isnan(result):
        result = 0.
    return result

def Y_sz(x):
    result = G_bb(x)*(x*(np.exp(x)+1.)/(np.exp(x)-1.)-4.)
    if math.isnan(result):
        result = 0.
    return result

def M(x):
    result = (x/beta_mu-1.)*np.exp(x)/(np.exp(x)-1.)**2.
    if math.isnan(result):
        result = 0.
    return result


def DI_normalization_in_MJy_per_sr(x,cosmo):
    nu_in_Hz = nu_in_GHz_of_x(x,cosmo)*1e9
    norm_in_Jy_per_sr = 2.*hplanck*nu_in_Hz**3./clight**2.*1e26
    return norm_in_Jy_per_sr*1e6



def nu_in_GHz_of_x(x,cosmo):
    return kb*cosmo.T_cmb*x/hplanck*1e-9
