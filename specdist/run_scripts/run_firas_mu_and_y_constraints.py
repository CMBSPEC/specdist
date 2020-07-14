import specdist as sd
from scipy.optimize import curve_fit
import numpy as np


T_pivot = 2.725
print('T0 = %.4e'%T_pivot)
firas = sd.firas()
sd.GetFirasXTatTpivot(firas,T_pivot)
firas.firas_x
# define function to fit firas monopole in MJy/sr
# the parameters are T, G0, mu, y
def cmb_monopole_linear_fit(x_array,T_cmb,G0,mu,y):
    delta_T = T_cmb - T_pivot
    #DI_ct = []
    monopole = []
    galaxy = []
    y_dist = []

    for i in range(len(firas.firas_x)):
        #yi = self.f(self.x_data[i])
        #DI_ct.append(finj*yi)
        nu = firas.firas_nu[i]
        monopole.append(sd.B_nu_of_T(nu,T_pivot)+delta_T*sd.dB_nu_dT_at_T(nu,T_pivot))
        galaxy.append(G0*firas.firas_Gnu[i])
        y_dist.append(y*2.*sd.hplanck*nu**3./sd.clight**2.*1e20*sd.Y_sz(firas.firas_x[i]))
        #galaxy.append(G0*nu**2.*sd.B_nu_of_T(nu,9.))

    #DI_ct = np.asarray(DI_ct)
    galaxy = np.asarray(galaxy)
    monopole = np.asarray(monopole)
    #mu_dist = sd.GetMuSpecDistAtTandX_chluba(mu,T_pivot,firas.firas_x)
    mu_dist = sd.GetMuSpecDistAtTandX(mu,T_pivot,firas.firas_x)
    #y_dist = sd.GetYSpecDistAtTandX(y,sd.firas_T0_bf,sd.firas.firas_x)
    y_dist = np.asarray(y_dist)
    return galaxy + monopole + y_dist + mu_dist


def cmb_monopole_nonlinear_fit(x_array,T_cmb,G0,mu,y):
    delta_T = 0.
    #DI_ct = []
    monopole = []
    galaxy = []
    y_dist = []

    for i in range(len(firas.firas_x)):
        #yi = self.f(self.x_data[i])
        #DI_ct.append(finj*yi)
        nu = firas.firas_nu[i]
        monopole.append(sd.B_nu_of_T(nu,T_cmb)+delta_T*sd.dB_nu_dT_at_T(nu,T_pivot))
        galaxy.append(G0*firas.firas_Gnu[i])
        y_dist.append(y*2.*sd.hplanck*nu**3./sd.clight**2.*1e20*sd.Y_sz(firas.firas_x[i]))
        #galaxy.append(G0*nu**2.*sd.B_nu_of_T(nu,9.))

    #DI_ct = np.asarray(DI_ct)
    galaxy = np.asarray(galaxy)
    monopole = np.asarray(monopole)
    #mu_dist = sd.GetMuSpecDistAtTandX_chluba(mu,T_pivot,firas.firas_x)
    mu_dist = sd.GetMuSpecDistAtTandX(mu,T_pivot,firas.firas_x)
    #y_dist = sd.GetYSpecDistAtTandX(y,sd.firas_T0_bf,sd.firas.firas_x)
    y_dist = np.asarray(y_dist)
    return galaxy + monopole + y_dist + mu_dist




def cmb_monopole_linear_fit_gal(x_array,T_cmb,G0,Tgal,mu,y):
    delta_T = T_cmb - T_pivot
    #DI_ct = []
    monopole = []
    galaxy = []
    y_dist = []

    for i in range(len(firas.firas_x)):
        #yi = self.f(self.x_data[i])
        #DI_ct.append(finj*yi)
        nu = firas.firas_nu[i]
        monopole.append(sd.B_nu_of_T(nu,T_pivot)+delta_T*sd.dB_nu_dT_at_T(nu,T_pivot))
        #galaxy.append(G0*firas.firas_Gnu[i])
        y_dist.append(y*2.*sd.hplanck*nu**3./sd.clight**2.*1e20*sd.Y_sz(firas.firas_x[i]))
        galaxy.append(G0*nu**2.*sd.B_nu_of_T(nu,Tgal/9.))

    #DI_ct = np.asarray(DI_ct)
    galaxy = np.asarray(galaxy)
    monopole = np.asarray(monopole)
    mu_dist = sd.GetMuSpecDistAtTandX_chluba(mu,sd.firas_T0_bf,sd.firas.firas_x)
    #mu_dist = sd.GetMuSpecDistAtTandX(mu,T_pivot,firas.firas_x)
    #y_dist = sd.GetYSpecDistAtTandX(y,sd.firas_T0_bf,sd.firas.firas_x)
    y_dist = np.asarray(y_dist)


    return galaxy + monopole + y_dist + mu_dist

X = firas.firas_x

# print('T_cmb -- G0 -- mu -- y')
# popt, pcov = curve_fit(cmb_monopole_linear_fit, X, sd.firas.firas_I0,sigma=sd.firas.firas_covmat,absolute_sigma=True)
# perr = np.sqrt(np.diag(pcov))
# print(popt)
# print(perr)


# print('T_cmb -- G0')
# popt,pcov = curve_fit(lambda x, Tcmb, G0, Tgal: cmb_monopole_linear_fit_gal(x, Tcmb, G0, Tgal, 0, 0),X, firas.firas_I0,sigma=firas.firas_covmat,absolute_sigma=True)
# perr = np.sqrt(np.diag(pcov))
# print(popt)
# print(perr)
# #rmu = 2.*np.sqrt(perr[2]**2.+firas.firas_mu_1996_systematic_stddev**2.)
# print("Including sys. |mu| < %.4e (95CL)" % rmu)

# bf_model = cmb_monopole_linear_fit_gal(X,popt[0],popt[1],popt[2],0.,0.)
# D = firas.firas_I0-bf_model
# chi2 = np.dot(D, np.dot(np.linalg.inv(firas.firas_covmat), D))
# print(chi2)




print('T_cmb -- G0 ')
popt,pcov = curve_fit(lambda x, Tcmb, G0: cmb_monopole_linear_fit(x, Tcmb, G0, 0, 0),X, firas.firas_I0,sigma=firas.firas_covmat,absolute_sigma=True)
perr = np.sqrt(np.diag(pcov))
print(popt)
print(perr)
# rmu = 2.*np.sqrt(perr[2]**2.+firas.firas_mu_1996_systematic_stddev**2.)
# print("Including sys. |mu| < %.4e (95CL)" % rmu)

bf_model = cmb_monopole_linear_fit(X,popt[0],popt[1],0.,0.)
D = firas.firas_I0-bf_model
chi2 = np.dot(D, np.dot(np.linalg.inv(firas.firas_covmat), D))
print('chi2/dof = %.5f'%(chi2/41.))


print('T_cmb -- mu ')
popt,pcov = curve_fit(lambda x, Tcmb, mu: cmb_monopole_linear_fit(x, Tcmb, -1e-3, mu, 0),X, firas.firas_I0,sigma=firas.firas_covmat,absolute_sigma=True)
perr = np.sqrt(np.diag(pcov))
print(popt)
print(perr)
rmu = 2.*np.sqrt(perr[1]**2.+firas.firas_mu_1996_systematic_stddev**2.)
print("Including sys. |mu| < %.4e (95CL)" % rmu)

bf_model = cmb_monopole_linear_fit(X,popt[0],-1e-3,popt[1],0.)
D = firas.firas_I0-bf_model
chi2 = np.dot(D, np.dot(np.linalg.inv(firas.firas_covmat), D))
print('chi2/dof = %.5f'%(chi2/41.))



print('T_cmb -- y ')
popt,pcov = curve_fit(lambda x, Tcmb, y: cmb_monopole_linear_fit(x, Tcmb, -1e-3, 0, y),X, firas.firas_I0,sigma=firas.firas_covmat,absolute_sigma=True)
perr = np.sqrt(np.diag(pcov))
print(popt)
print(perr)
rmu = 2.*np.sqrt(perr[1]**2.+firas.firas_mu_1996_systematic_stddev**2.)
print("Including sys. |y| < %.4e (95CL)" % rmu)

bf_model = cmb_monopole_linear_fit(X,popt[0],-1e-3,0.,popt[1])
D = firas.firas_I0-bf_model
chi2 = np.dot(D, np.dot(np.linalg.inv(firas.firas_covmat), D))
print('chi2/dof = %.5f'%(chi2/41.))




print('T_cmb -- G0 -- mu')
popt,pcov = curve_fit(lambda x, Tcmb, G0, mu: cmb_monopole_linear_fit(x, Tcmb, G0, mu, 0),X, firas.firas_I0,sigma=firas.firas_covmat,absolute_sigma=True)
perr = np.sqrt(np.diag(pcov))
print(popt)
print(perr)
rmu = 2.*np.sqrt(perr[2]**2.+firas.firas_mu_1996_systematic_stddev**2.)
print("Including sys. |mu| < %.4e (95CL)" % rmu)

bf_model = cmb_monopole_linear_fit(X,popt[0],popt[1],popt[2],0.)
D = firas.firas_I0-bf_model
chi2 = np.dot(D, np.dot(np.linalg.inv(firas.firas_covmat), D))
print('chi2/dof = %.5f'%(chi2/40.))


print('T_cmb -- G0 -- y')
popt,pcov = curve_fit(lambda x, Tcmb, G0, y: cmb_monopole_linear_fit(x, Tcmb, G0, 0, y),X, firas.firas_I0,sigma=firas.firas_covmat,absolute_sigma=True)
perr = np.sqrt(np.diag(pcov))
print(popt)
print(perr)
rmu = 2.*np.sqrt(perr[2]**2.+firas.firas_y_1996_systematic_stddev**2.)
print("Including sys. |y| < %.4e (95CL)" % rmu)

bf_model = cmb_monopole_linear_fit(X,popt[0],popt[1],0.,popt[2])
D = firas.firas_I0-bf_model
chi2 = np.dot(D, np.dot(np.linalg.inv(firas.firas_covmat), D))
print('chi2/dof = %.5f'%(chi2/40.))



print('NL --- T_cmb -- G0 -- mu')
popt,pcov = curve_fit(lambda x, Tcmb, G0, mu: cmb_monopole_nonlinear_fit(x, Tcmb, G0, mu, 0),X, firas.firas_I0,sigma=firas.firas_covmat,absolute_sigma=True)
perr = np.sqrt(np.diag(pcov))
print(popt)
print(perr)
rmu = 2.*np.sqrt(perr[2]**2.+firas.firas_mu_1996_systematic_stddev**2.)
print("Including sys. |mu| < %.4e (95CL)" % rmu)

bf_model = cmb_monopole_linear_fit(X,popt[0],popt[1],popt[2],0.)
D = firas.firas_I0-bf_model
chi2 = np.dot(D, np.dot(np.linalg.inv(firas.firas_covmat), D))
print('chi2/dof = %.5f'%(chi2/40.))


print('NL --- T_cmb -- G0 -- y')
popt,pcov = curve_fit(lambda x, Tcmb, G0, y: cmb_monopole_nonlinear_fit(x, Tcmb, G0, 0, y),X, firas.firas_I0,sigma=firas.firas_covmat,absolute_sigma=True)
perr = np.sqrt(np.diag(pcov))
print(popt)
print(perr)
rmu = 2.*np.sqrt(perr[2]**2.+firas.firas_y_1996_systematic_stddev**2.)
print("Including sys. |y| < %.4e (95CL)" % rmu)

bf_model = cmb_monopole_linear_fit(X,popt[0],popt[1],0.,popt[2])
D = firas.firas_I0-bf_model
chi2 = np.dot(D, np.dot(np.linalg.inv(firas.firas_covmat), D))
print('chi2/dof = %.5f'%(chi2/40.))


print('T_cmb -- G0 -- mu -- y')
popt,pcov = curve_fit(lambda x, Tcmb, G0, mu, y: cmb_monopole_linear_fit(x, Tcmb, G0, mu, y),X, firas.firas_I0,sigma=firas.firas_covmat,absolute_sigma=True)
perr = np.sqrt(np.diag(pcov))
print(popt)
print(perr)
rmu = 2.*np.sqrt(perr[2]**2.+firas.firas_mu_1996_systematic_stddev**2.)
ry = 2.*np.sqrt(perr[3]**2.+firas.firas_y_1996_systematic_stddev**2.)
print("Including sys. |mu| < %.4e (95CL)" % rmu)
print("Including sys. |y| < %.4e (95CL)" % ry)

bf_model = cmb_monopole_linear_fit(X,popt[0],popt[1],0.,popt[2])
D = firas.firas_I0-bf_model
chi2 = np.dot(D, np.dot(np.linalg.inv(firas.firas_covmat), D))
print('chi2/dof = %.5f'%(chi2/39.))


# print('T_cmb -- G0 -- y')
# popt,pcov = curve_fit(lambda x, Tcmb, G0, y: cmb_monopole_linear_fit(x, Tcmb, G0, 0, y),X, sd.firas.firas_I0,sigma=sd.firas.firas_covmat,absolute_sigma=True)
# perr = np.sqrt(np.diag(pcov))
# print(popt)
# print(perr)
# rmu = 2.*np.sqrt(perr[2]**2.+sd.firas.firas_y_1996_systematic_stddev**2.)
# print("Including sys. |y| < %.4e (95CL)" % rmu)
