import specdist as sd
from scipy.optimize import curve_fit
import numpy as np


#print(sd.firas.firas_Q)
#print(sd.firas.firas_covmat)
#print(sd.firas.nu_minus_nuprime)
#exit(0)



# define function to fit firas monopole in MJy/sr
# the parameters are T, G0, mu, y
def cmb_monopole_linear_fit(x_array,T_cmb,G0,mu,y):
    delta_T = T_cmb - sd.firas_T0_bf
    #DI_ct = []
    monopole = []
    galaxy = []
    y_dist = []

    for i in range(len(sd.firas.firas_x)):
        #yi = self.f(self.x_data[i])
        #DI_ct.append(finj*yi)
        nu = sd.firas.firas_nu[i]
        monopole.append(sd.B_nu_of_T(nu,sd.firas_T0_bf)+delta_T*sd.dB_nu_dT_at_T(nu,sd.firas_T0_bf))
        galaxy.append(G0*sd.firas.firas_Gnu[i])
        y_dist.append(y*2.*sd.hplanck*nu**3./sd.clight**2.*1e20*sd.Y_sz(sd.firas.firas_x[i]))
        #galaxy.append(G0*nu**2.*sd.B_nu_of_T(nu,9.))

    #DI_ct = np.asarray(DI_ct)
    galaxy = np.asarray(galaxy)
    monopole = np.asarray(monopole)
    #mu_dist = sd.GetMuSpecDistAtTandX_chluba(mu,sd.firas_T0_bf,sd.firas.firas_x)
    mu_dist = sd.GetMuSpecDistAtTandX(mu,sd.firas_T0_bf,sd.firas.firas_x)
    #y_dist = sd.GetYSpecDistAtTandX(y,sd.firas_T0_bf,sd.firas.firas_x)
    y_dist = np.asarray(y_dist)


    return galaxy + monopole + y_dist + mu_dist

X = sd.firas.firas_x

print('T_cmb -- G0 -- mu -- y')
popt, pcov = curve_fit(cmb_monopole_linear_fit, X, sd.firas.firas_I0,sigma=sd.firas.firas_covmat,absolute_sigma=True)
perr = np.sqrt(np.diag(pcov))
print(popt)
print(perr)

print('T_cmb -- G0 -- mu')
popt,pcov = curve_fit(lambda x, Tcmb, G0, mu: cmb_monopole_linear_fit(x, Tcmb, G0, mu, 0),X, sd.firas.firas_I0,sigma=sd.firas.firas_covmat,absolute_sigma=True)
perr = np.sqrt(np.diag(pcov))
print(popt)
print(perr)
rmu = 2.*np.sqrt(perr[2]**2.+sd.firas.firas_mu_1996_systematic_stddev**2.)
print("Including sys. |mu| < %.4e (95CL)" % rmu)

print('T_cmb -- G0 -- y')
popt,pcov = curve_fit(lambda x, Tcmb, G0, y: cmb_monopole_linear_fit(x, Tcmb, G0, 0, y),X, sd.firas.firas_I0,sigma=sd.firas.firas_covmat,absolute_sigma=True)
perr = np.sqrt(np.diag(pcov))
print(popt)
print(perr)
rmu = 2.*np.sqrt(perr[2]**2.+sd.firas.firas_y_1996_systematic_stddev**2.)
print("Including sys. |y| < %.4e (95CL)" % rmu)

photon_injection ='yes'
photon_injection_case = 'bare'

if photon_injection == 'yes':
    sd_lib = sd.specdist_ct_spectra_lib()
    sd.load_ct_spectra_lib(photon_injection_case,sd_lib)

def pi_cmb_monopole_linear_fit(x_array,T_cmb,G0,finj_over_fct,gamma_asked,xinj_asked):
    delta_T = T_cmb - sd.firas_T0_bf
    #DI_ct = []
    monopole = []
    galaxy = []


    S = sd.GetSpectra(gamma_asked,xinj_asked,sd.firas.firas_x,sd_lib)
    pispec = finj_over_fct*np.asarray(S["DI"])*1e-6

    for i in range(len(sd.firas.firas_x)):
        #yi = self.f(self.x_data[i])
        #DI_ct.append(finj*yi)
        nu = sd.firas.firas_nu[i]
        monopole.append(sd.B_nu_of_T(nu,sd.firas_T0_bf)+delta_T*sd.dB_nu_dT_at_T(nu,sd.firas_T0_bf))
        galaxy.append(G0*sd.firas.firas_Gnu[i])
        #galaxy.append(G0*nu**2.*sd.B_nu_of_T(nu,9.))


    #DI_ct = np.asarray(DI_ct)
    galaxy = np.asarray(galaxy)
    monopole = np.asarray(monopole)
    return galaxy + monopole + pispec

print('T_cmb -- G0 -- finj_over_fct')
for G_X in np.logspace(-17,-8,10):
    for xinj in np.logspace(-6,6,10):
        popt,pcov = curve_fit(lambda x, Tcmb, G0, finj_over_fct: pi_cmb_monopole_linear_fit(x, Tcmb, G0, finj_over_fct,G_X,xinj),X, sd.firas.firas_I0,sigma=sd.firas.firas_covmat,absolute_sigma=True)
        perr = np.sqrt(np.diag(pcov))
        #print(popt)
        print(G_X,xinj,perr[2])
# rmu = 2.*np.sqrt(perr[2]**2.+sd.firas.firas_y_1996_systematic_stddev**2.)
# print("Including sys. |y| < %.4e (95CL)" % rmu)
