from .utils import *
from pkg_resources import resource_filename
import numpy as np
import re
from scipy.interpolate import interp1d
import os





def rho_gamma_in_GeV_per_cm3(T):
    return 8.*np.pi**5.*(kb*T)**4/15./clight**3./hplanck**3./GeV_over_kg*1e-6/clight**2

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
