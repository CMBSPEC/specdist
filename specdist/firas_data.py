from .utils import *

#D = np.loadtxt(self.firas_data_file)

class firas:
    def __init__(self):
        self.firas_data_file = resource_filename("specdist","data/firas_monopole_spec_v1.txt")


        self.firas_x = np.array([], 'float64')
        self.firas_nu = np.array([], 'float64')
        self.firas_I0 = np.array([], 'float64')
        self.firas_residual = np.array([], 'float64')
        self.firas_sigma = np.array([], 'float64')
        self.firas_Gnu = np.array([], 'float64')

        for line in open(self.firas_data_file, 'r'):
            if (line.find('#') == -1):
                self.firas_nu = np.append(self.firas_nu, float(line.split()[0]))
                self.firas_I0 = np.append(self.firas_I0, float(line.split()[1]))
                self.firas_residual = np.append(self.firas_residual, float(line.split()[2]))
                self.firas_sigma = np.append(self.firas_sigma, float(line.split()[3])*1e-3)
                self.firas_Gnu = np.append(self.firas_Gnu, float(line.split()[4])*1e-3)
        self.firas_nu_cm = self.firas_nu
        self.firas_nu = np.asarray(self.firas_nu)*clight*1.e2 # in s^-1



        self.firas_covmat = np.zeros((len(self.firas_nu),len(self.firas_nu)))
        nu_minus_nuprime = self.firas_nu_cm - self.firas_nu_cm[0]
        self.firas_Q = np.array([1.000, 0.176,-0.203, 0.145, 0.077,-0.005,-0.022, 0.032, 0.053, 0.025,-0.003, 0.007, 0.029, 0.029, 0.003,-0.002, 0.016, 0.020, 0.011, 0.002, 0.007, 0.011, 0.009, 0.003,-0.004,-0.001, 0.003, 0.003,-0.001,-0.003, 0.000, 0.003, 0.009, 0.015, 0.008, 0.003,-0.002, 0.000,-0.006,-0.006, 0.000, 0.002, 0.008])

        #print(nu_minus_nuprime)

        for inu in range(len(self.firas_nu)):
            for inuprime in range(len(self.firas_nu)):
                nu  = self.firas_nu_cm[inu]
                nuprime  = self.firas_nu_cm[inuprime]
                r = round(np.abs(nu-nuprime),4)
                iq = (np.abs(nu_minus_nuprime - r)).argmin()
                #iq = np.where(np.roll(nu_minus_nuprime,inu) == r)
                #print(inu,inuprime,iq)
                self.firas_covmat[inu][inuprime] = self.firas_sigma[inu]*self.firas_sigma[inuprime]*self.firas_Q[iq]


        self.firas_nu_min_in_GHz = np.min(self.firas_nu)/1e9
        self.firas_nu_max_in_GHz = np.max(self.firas_nu)/1e9



        self.firas_mu_1996_95_cl = 9.e-5
        self.firas_y_1996_95_cl = 15.e-6

        self.firas_mu_1996_systematic_stddev = 1e-5
        self.firas_y_1996_systematic_stddev = 4e-6

        self.firas_Drho_rho_1996_95_cl = 6e-5 #approx 4*firas_y_1996=6e-05, firas_mu_1996/1.401

        GetFirasXTatTpivot(self,firas_T0_bf)

def GetFirasXTatTpivot(firas,Tpivot):
    firas.firas_x = hplanck*firas.firas_nu/kb/Tpivot
    firas.firas_x_min = np.min(firas.firas_x)
    firas.firas_x_max = np.max(firas.firas_x)
