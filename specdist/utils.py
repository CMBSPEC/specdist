import numpy as np
from datetime import datetime
import multiprocessing
import functools
import re
from pkg_resources import resource_filename
import os
from scipy import optimize
from scipy.integrate import quad
from scipy.interpolate import interp1d
import math


def find_nearests(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    if array[idx]>value:
        idxp =idx
        idxm = idx-1
    else:
        idxm = idx
        idxp = idx+1
    if idxp == len(array):
        idxm -= 1
        idxp -= 1
    return (idxm,idxp)

def scientific_notation(p_value,digit=2):
    str_xinj_asked = str("%.3e"%p_value)
    text_gamma_str1 = ''
    if p_value>1.:
        num = float(str_xinj_asked.split('e+')[0])
        exp = int(str_xinj_asked.split('e+')[1])
        if (round(num,digit)==10.):
            num = 1.
            exp = exp + 1
        if digit == 1:
            text_gamma_str1 = r'$%.1f \times 10^{%d}$'% (num,exp)
        if digit == 2:
            text_gamma_str1 = r'$%.2f \times 10^{%d}$'% (num,exp)
        if digit == 0:
            text_gamma_str1 = r'$%.0f \times 10^{%d}$'% (num,exp)
        if num == 1.:
            text_gamma_str1 = r'$10^{%d}$'% (exp)
    if p_value<1.:
        num = float(str_xinj_asked.split('e-')[0])
        exp = int(str_xinj_asked.split('e-')[1])
        if (round(num,digit)==10.):
            num = 1.
            exp = exp - 1
        if digit == 1:
            text_gamma_str1 = r'$%.1f \times 10^{-%d}$'% (num,exp)
        if digit == 2:
            text_gamma_str1 = r'$%.2f \times 10^{-%d}$'% (num,exp)
        if digit == 0:
            text_gamma_str1 = r'$%.0f \times 10^{-%d}$'% (num,exp)
        if num == 1.:
            text_gamma_str1 = r'$10^{-%d}$'% (exp)
    if p_value==1.:
        text_gamma_str1 = r'$1$'
    return text_gamma_str1





#1 GeV/c2 = 1.78266192×10−27 kg
GeV_over_kg = 1.78266192e-27

#1 km/Mpc
km_over_Mpc = 3.24077929e-20

# 1 GHz/eV
GHz_over_eV =  4.1356655385e-6



def nu_in_GHz_of_x(x,cosmo):
    return kb*cosmo.T_cmb*x/hplanck*1e-9

def x_of_nu_in_GHz(nu,cosmo):
    return hplanck*nu*1e9/kb/cosmo.T_cmb

def x_of_hnu_in_eV(hnu,cosmo):
    nu_in_GHz = hnu*1./GHz_over_eV
    return x_of_nu_in_GHz(nu_in_GHz,cosmo)

kb = 1.38064852e-23 #m2 kg s-2 K-1
clight = 299792458. #m/s
hplanck=6.62607004e-34 #m2 kg / s
firas_T0 = 2.728 #pivot temperature used in the Max Lkl Analysis
firas_T0_bf = 2.725 #best-fitting temperature
G_newton = 6.674e-11
rho_crit_over_h2_in_GeV_per_cm3 = 1.0537e-5


nu_21_cm_in_GHz =  1./21.1*clight*1.e2/1.e9
x_21_cm = hplanck*nu_21_cm_in_GHz/kb/firas_T0_bf*1.e9

kappa_c = 2.1419 # 4M_2-3M_c see below eq. 9b of https://arxiv.org/pdf/1506.06582.pdf

beta_mu = 2.1923

G1 = np.pi**2./6
G2 = 2.4041
G3 = np.pi**4/15.
a_rho = G2/G3
alpha_mu = 2.*G1/3./G2 # = 1/beta_mu = π^2/18ζ(3) see eq. 4.15 CUSO lectures.

z_mu_era = 3e5
z_y_era = 5e4
z_reio_min = 6
z_reio_max = 25
z_recombination_min = 800
z_recombination_max = 1500 
