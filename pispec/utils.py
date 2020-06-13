import numpy as np
from datetime import datetime
import multiprocessing
import functools
import re
from pkg_resources import resource_filename

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

kb = 1.38064852e-23 #m2 kg s-2 K-1
clight = 299792458. #m/s
hplanck=6.62607004e-34 #m2 kg / s
firas_T0 = 2.728 #pivot temperature used in the Max Lkl Analysis
firas_T0_bf = 2.725 #best-fitting temperature
rho_crit_in_h2_Gev_per_cm3 = 1.0537e-5
nu_21_cm_in_GHz =  1./21.1*clight*1.e2/1.e9
x_21_cm = hplanck*nu_21_cm_in_GHz/kb/firas_T0_bf*1.e9
