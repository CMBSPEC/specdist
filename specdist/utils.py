import numpy as np
from datetime import datetime
import multiprocessing
import functools
import re
from pkg_resources import resource_filename
import os
from scipy import optimize
from scipy.integrate import quad


def find_nearests(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    if array[idx]>value:
        idxp =idx
        idxm = idx-1
    else:
        idxm = idx
        idxp = idx+1
    return (idxm,idxp)


#1 GeV/c2 = 1.78266192×10−27 kg
GeV_over_kg = 1.78266192e-27

#1 km/Mpc
km_over_Mpc = 3.24077929e-20

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


G2 = 2.4041
G3 = np.pi**4/15.
a_rho = G2/G3
