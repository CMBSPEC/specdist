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
            dist.append(mu*(-dS_dMU(nu,T)*xp*(0.4561-1./xp)))
    except:
        nu = kb*T/hplanck*x
        dist.append(mu*(-dS_dMU(nu,T)*x*(0.4561-1./x)))
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
    return x*np.exp(x)/(np.exp(x)-1.)**2.

def Y_sz(x):
    return G_bb(x)*(x*(np.exp(x)+1.)/(np.exp(x)-1.)-4.)
