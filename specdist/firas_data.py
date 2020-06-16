from .utils import *

#D = np.loadtxt(self.firas_data_file)

class firas:
    firas_data_file = resource_filename("pispec","data/firas_monopole_spec_v1.txt")


    firas_x = np.array([], 'float64')
    firas_nu = np.array([], 'float64')
    firas_I0 = np.array([], 'float64')
    firas_residual = np.array([], 'float64')
    firas_sigma = np.array([], 'float64')
    firas_Gnu = np.array([], 'float64')

    for line in open(firas_data_file, 'r'):
        if (line.find('#') == -1):
            firas_nu = np.append(firas_nu, float(line.split()[0]))
            firas_I0 = np.append(firas_I0, float(line.split()[1]))
            firas_residual = np.append(firas_residual, float(line.split()[2]))
            firas_sigma = np.append(firas_sigma, float(line.split()[3])*1e-3)
            firas_Gnu = np.append(firas_Gnu, float(line.split()[4])*1e-3)
    firas_nu = np.asarray(firas_nu)*clight*1.e2 # in s^-1
    firas_x = hplanck*firas_nu/kb/firas_T0_bf


    firas_nu_min_in_GHz = np.min(firas_nu)/1e9
    firas_nu_max_in_GHz = np.max(firas_nu)/1e9

    firas_x_min = np.min(firas_x)
    firas_x_max = np.max(firas_x)

    firas_mu_1996_95_cl = 9.e-5
    firas_y_1996_95_cl = 15.e-6

    firas_Drho_rho_1996_95_cl = 6e-5 #approx 4*firas_y_1996=6e-05, firas_mu_1996/1.401
