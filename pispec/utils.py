import numpy as np


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



kb = 1.38064852e-23
clight = 299792458.
hplanck=6.62607004e-34
firas_T0 = 2.728 #pivot temperature used in the Max Lkl Analysis
firas_T0_bf = 2.725 #best-fitting temperature

nu_21_cm_in_GHz =  1./21.1*clight*1.e2/1.e9
x_21_cm = hplanck*nu_21_cm_in_GHz/kb/firas_T0_bf*1.e9
