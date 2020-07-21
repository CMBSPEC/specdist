from .utils import *
from .config import *
# from pkg_resources import resource_filename
# import numpy as np
# import re
# from scipy.interpolate import interp1d
# import os


class specdist_ct_spectra_lib:
    def __init__(self):
        self.case_id = ""

        # this_dir, this_filename = os.path.split(__file__)
        # DATA_PATH = os.path.join(this_dir, "data", "data.txt")
        # path_to_spectra = resource_filename("specdist","data/ct_database/"+case_id)
        # print(path_to_spectra)
        self.path_to_spectra = ""

        self.Gamma_inj_min = 1e-17
        self.Gamma_inj_max = 1e-8
        self.N_Gamma_inj = 50
        self.Gamma_values = np.logspace(np.log10(self.Gamma_inj_min),np.log10(self.Gamma_inj_max),self.N_Gamma_inj)

        self.x_inj_min = 1e-6
        self.x_inj_max = 1e6
        self.N_x_inj = 200
        self.x_inj_values = np.logspace(np.log10(self.x_inj_min),np.log10(self.x_inj_max),self.N_x_inj)


        self.X_2d = []
        self.DI_2d = []
        self.finj_2d = []

        self.Xe_values_2d = []
        self.Xe_values_no_inj_2d = []
        self.Xe_redshifts_2d = []
        self.Xe_redshifts_no_inj_2d = []

        self.DXe_Xe_2d = []
        self.DXe_Xe_redshifts_2d = []
        self.Xe_PCA_EigenModes = {}
        self.Xe_PCA_EigenModes['E1'] = {}
        self.Xe_PCA_EigenModes['E1']['z'] = []
        self.Xe_PCA_EigenModes['E1']['values'] = []
        self.Xe_PCA_EigenModes['E2'] = {}
        self.Xe_PCA_EigenModes['E2']['z'] = []
        self.Xe_PCA_EigenModes['E2']['values'] = []
        self.Xe_PCA_EigenModes['E3'] = {}
        self.Xe_PCA_EigenModes['E3']['z'] = []
        self.Xe_PCA_EigenModes['E3']['values'] = []



def load_ct_spectra_lib(case,specdist_ct_spectra_lib):
    if case == 'bare':
        specdist_ct_spectra_lib.case_id = "case_1_040520"
    elif case == 'lyc':
        specdist_ct_spectra_lib.case_id = "case_lyman_090620"
    elif case == 'lyc_reio':
        specdist_ct_spectra_lib.case_id = "case_lyman_reio_180620"
    elif case == 'raw_lyc_reio':
        specdist_ct_spectra_lib.case_id = "case_raw_lyman_reio_180620"
    elif case == 'mu_fit':
        specdist_ct_spectra_lib.case_id = "case_for_mu_fit"
        specdist_ct_spectra_lib.Gamma_inj_min = 1e-9
        specdist_ct_spectra_lib.Gamma_inj_max = 1e-6
        specdist_ct_spectra_lib.N_Gamma_inj = 4
        specdist_ct_spectra_lib.Gamma_values = np.logspace(np.log10(specdist_ct_spectra_lib.Gamma_inj_min),np.log10(specdist_ct_spectra_lib.Gamma_inj_max),specdist_ct_spectra_lib.N_Gamma_inj)

        specdist_ct_spectra_lib.x_inj_min = 1e1
        specdist_ct_spectra_lib.x_inj_max = 1e8
        specdist_ct_spectra_lib.N_x_inj = 50
        specdist_ct_spectra_lib.x_inj_values = np.logspace(np.log10(specdist_ct_spectra_lib.x_inj_min),np.log10(specdist_ct_spectra_lib.x_inj_max),specdist_ct_spectra_lib.N_x_inj)

    elif case == 'xe_history_200720':
        specdist_ct_spectra_lib.case_id = "case_" + case
        specdist_ct_spectra_lib.Gamma_inj_min = 1e-17
        specdist_ct_spectra_lib.Gamma_inj_max = 1e-12
        specdist_ct_spectra_lib.N_Gamma_inj = 10
        specdist_ct_spectra_lib.Gamma_values = np.logspace(np.log10(specdist_ct_spectra_lib.Gamma_inj_min),np.log10(specdist_ct_spectra_lib.Gamma_inj_max),specdist_ct_spectra_lib.N_Gamma_inj)

        specdist_ct_spectra_lib.x_inj_min = 1e-8
        specdist_ct_spectra_lib.x_inj_max = 1e7
        specdist_ct_spectra_lib.N_x_inj = 200
        specdist_ct_spectra_lib.x_inj_values = np.logspace(np.log10(specdist_ct_spectra_lib.x_inj_min),np.log10(specdist_ct_spectra_lib.x_inj_max),specdist_ct_spectra_lib.N_x_inj)
        E1 = np.loadtxt(path_to_ct_database+'../PCA_modes/Modes/mode_N121_so_planck_1.dat')
        E2 = np.loadtxt(path_to_ct_database+'../PCA_modes/Modes/mode_N121_so_planck_2.dat')
        E3 = np.loadtxt(path_to_ct_database+'../PCA_modes/Modes/mode_N121_so_planck_3.dat')
        specdist_ct_spectra_lib.Xe_PCA_EigenModes['E1']['z'] = E1[:,0]
        specdist_ct_spectra_lib.Xe_PCA_EigenModes['E2']['z'] = E2[:,0]
        specdist_ct_spectra_lib.Xe_PCA_EigenModes['E3']['z'] = E3[:,0]
        specdist_ct_spectra_lib.Xe_PCA_EigenModes['E1']['values'] = E1[:,1]
        specdist_ct_spectra_lib.Xe_PCA_EigenModes['E2']['values'] = E2[:,1]
        specdist_ct_spectra_lib.Xe_PCA_EigenModes['E3']['values'] = E3[:,1]

    elif case == 'xe_history_200720_finj_fisher':
        specdist_ct_spectra_lib.case_id = "case_" + case
        specdist_ct_spectra_lib.Gamma_inj_min = 1e-17
        specdist_ct_spectra_lib.Gamma_inj_max = 1e-12
        specdist_ct_spectra_lib.N_Gamma_inj = 1
        specdist_ct_spectra_lib.Gamma_values = np.logspace(np.log10(specdist_ct_spectra_lib.Gamma_inj_min),np.log10(specdist_ct_spectra_lib.Gamma_inj_max),specdist_ct_spectra_lib.N_Gamma_inj)

        specdist_ct_spectra_lib.x_inj_min = 1e-8
        specdist_ct_spectra_lib.x_inj_max = 1e7
        specdist_ct_spectra_lib.N_x_inj = 200
        specdist_ct_spectra_lib.x_inj_values = np.logspace(np.log10(specdist_ct_spectra_lib.x_inj_min),np.log10(specdist_ct_spectra_lib.x_inj_max),specdist_ct_spectra_lib.N_x_inj)
        E1 = np.loadtxt(path_to_ct_database+'../PCA_modes/Modes/mode_N121_so_planck_1.dat')
        E2 = np.loadtxt(path_to_ct_database+'../PCA_modes/Modes/mode_N121_so_planck_2.dat')
        E3 = np.loadtxt(path_to_ct_database+'../PCA_modes/Modes/mode_N121_so_planck_3.dat')
        specdist_ct_spectra_lib.Xe_PCA_EigenModes['E1']['z'] = E1[:,0]
        specdist_ct_spectra_lib.Xe_PCA_EigenModes['E2']['z'] = E2[:,0]
        specdist_ct_spectra_lib.Xe_PCA_EigenModes['E3']['z'] = E3[:,0]
        specdist_ct_spectra_lib.Xe_PCA_EigenModes['E1']['values'] = E1[:,1]
        specdist_ct_spectra_lib.Xe_PCA_EigenModes['E2']['values'] = E2[:,1]
        specdist_ct_spectra_lib.Xe_PCA_EigenModes['E3']['values'] = E3[:,1]

    elif 'xe_history' in case:
        specdist_ct_spectra_lib.case_id = "case_" + case
        specdist_ct_spectra_lib.Gamma_inj_min = 1e-17
        specdist_ct_spectra_lib.Gamma_inj_max = 1e-12
        specdist_ct_spectra_lib.N_Gamma_inj = 10
        specdist_ct_spectra_lib.Gamma_values = np.logspace(np.log10(specdist_ct_spectra_lib.Gamma_inj_min),np.log10(specdist_ct_spectra_lib.Gamma_inj_max),specdist_ct_spectra_lib.N_Gamma_inj)

        specdist_ct_spectra_lib.x_inj_min = 1e-8
        specdist_ct_spectra_lib.x_inj_max = 1e7
        specdist_ct_spectra_lib.N_x_inj = 64
        specdist_ct_spectra_lib.x_inj_values = np.logspace(np.log10(specdist_ct_spectra_lib.x_inj_min),np.log10(specdist_ct_spectra_lib.x_inj_max),specdist_ct_spectra_lib.N_x_inj)

    else:
        print('this case has not been computed. Computed cases are "lyc" or "bare".')
        return
    specdist_ct_spectra_lib.path_to_spectra  = path_to_ct_database + specdist_ct_spectra_lib.case_id
    if 'xe_history' in case :
        specdist_ct_spectra_lib.path_to_spectra  = path_to_ct_database + specdist_ct_spectra_lib.case_id + '/' + case

    if case == 'raw_lyc_reio':
        specdist_ct_spectra_lib.case_id = "case_lyman_reio_180620"
    if case == 'mu_fit':
        specdist_ct_spectra_lib.case_id = "for_mu_fit"



    for id_Gamma in range(specdist_ct_spectra_lib.N_Gamma_inj):
    # for id_Gamma in range(1):
        Gi = specdist_ct_spectra_lib.Gamma_values[id_Gamma]
        str_gamma = str("%.3e"%Gi)

        #read x array
        x_ct = []
        if case == 'mu_fit':
            filename = specdist_ct_spectra_lib.path_to_spectra + '/spectra_' + specdist_ct_spectra_lib.case_id + '_G_' + str_gamma +'/spectra_spectra_' + specdist_ct_spectra_lib.case_id + '_G_' + str_gamma  + '_x_ct.txt'
        elif 'xe_history' in case:
            filename = specdist_ct_spectra_lib.path_to_spectra + '_G_' + str_gamma + '/spectra_' + case + '_G_' + str_gamma + '_x_ct.txt'
        else:
            filename = specdist_ct_spectra_lib.path_to_spectra + '/spectra_' + specdist_ct_spectra_lib.case_id + '_G_' + str_gamma + '_x_ct.txt'
        with open(filename) as f:
            for line in f:
                ls = line.strip()
                if ls:
                    if "#" in ls:
                        continue
                    else:
                        x_ct_p = []
                        l = re.split('\t',ls)
                        l = [e for e in l if e]
                        #print(l)
                        for s in l:
                            x_cti = float(s)
                            x_ct_p.append(x_cti)
                        x_ct_p = np.asarray(x_ct_p)
                    x_ct.append(x_ct_p)


        #read DI array
        DI_ct = []
        if case == 'mu_fit':
            DI_ct_bare = []
            DI_ct_hubble = []
            filename = specdist_ct_spectra_lib.path_to_spectra + '/spectra_' + specdist_ct_spectra_lib.case_id + '_G_' + str_gamma +'/spectra_spectra_' + specdist_ct_spectra_lib.case_id + '_G_' + str_gamma  + '_DI_ct.txt'
            filename_hubble = specdist_ct_spectra_lib.path_to_spectra + '/spectra_' + specdist_ct_spectra_lib.case_id + '_hubble_G_' + str_gamma +'/spectra_spectra_' + specdist_ct_spectra_lib.case_id + '_hubble_G_' + str_gamma  + '_DI_ct.txt'
            with open(filename) as f:
                for line in f:
                    ls = line.strip()
                    if ls:
                        if "#" in ls:
                            continue
                        else:
                            DI_ct_p = []
                            l = re.split('\t',ls)
                            l = [e for e in l if e]
                            #print(l)
                            for s in l:
                                DI_cti = float(s)
                                DI_ct_p.append(DI_cti)
                            DI_ct_p = np.asarray(DI_ct_p)
                        DI_ct_bare.append(DI_ct_p)
            with open(filename_hubble) as f:
                for line in f:
                    ls = line.strip()
                    if ls:
                        if "#" in ls:
                            continue
                        else:
                            DI_ct_p = []
                            l = re.split('\t',ls)
                            l = [e for e in l if e]
                            #print(l)
                            for s in l:
                                DI_cti = float(s)
                                DI_ct_p.append(DI_cti)
                            DI_ct_p = np.asarray(DI_ct_p)
                        DI_ct_hubble.append(DI_ct_p)
            for (p_bare,p_hubble) in zip(DI_ct_bare,DI_ct_hubble):
                DI_ct.append(p_bare-p_hubble)
        elif 'xe_history' in case:
            DI_ct_bare = []
            DI_ct_hubble = []
            filename = specdist_ct_spectra_lib.path_to_spectra + '_G_' + str_gamma + '/spectra_' + case + '_G_' + str_gamma  + '_DI_ct.txt'
            filename_hubble = specdist_ct_spectra_lib.path_to_spectra + '_hubble_G_' + str_gamma + '/spectra_'  + case + '_hubble_G_' + str_gamma  + '_DI_ct.txt'
            with open(filename) as f:
                for line in f:
                    ls = line.strip()
                    if ls:
                        if "#" in ls:
                            continue
                        else:
                            DI_ct_p = []
                            l = re.split('\t',ls)
                            l = [e for e in l if e]
                            #print(l)
                            for s in l:
                                DI_cti = float(s)
                                DI_ct_p.append(DI_cti)
                            DI_ct_p = np.asarray(DI_ct_p)
                        DI_ct_bare.append(DI_ct_p)
            with open(filename_hubble) as f:
                for line in f:
                    ls = line.strip()
                    if ls:
                        if "#" in ls:
                            continue
                        else:
                            DI_ct_p = []
                            l = re.split('\t',ls)
                            l = [e for e in l if e]
                            #print(l)
                            for s in l:
                                DI_cti = float(s)
                                DI_ct_p.append(DI_cti)
                            DI_ct_p = np.asarray(DI_ct_p)
                        DI_ct_hubble.append(DI_ct_p)
            for (p_bare,p_hubble) in zip(DI_ct_bare,DI_ct_hubble):
                DI_ct.append(p_bare-p_hubble)
        else:
            filename = specdist_ct_spectra_lib.path_to_spectra + '/spectra_' + specdist_ct_spectra_lib.case_id + '_G_' + str_gamma + '_DI_ct.txt'
            with open(filename) as f:
                for line in f:
                    ls = line.strip()
                    if ls:
                        if "#" in ls:
                            continue
                        else:
                            DI_ct_p = []
                            l = re.split('\t',ls)
                            l = [e for e in l if e]
                            #print(l)
                            for s in l:
                                DI_cti = float(s)
                                DI_ct_p.append(DI_cti)
                            DI_ct_p = np.asarray(DI_ct_p)
                        DI_ct.append(DI_ct_p)



        if 'xe_history' in case:
            Xe_values_ct = []
            Xe_values_no_inj_ct = []
            Xe_redshifts_ct = []
            Xe_redshifts_no_inj_ct = []
            filename = specdist_ct_spectra_lib.path_to_spectra + '_G_' + str_gamma + '/spectra_' + case + '_G_' + str_gamma  + '_Xe_values_ct.txt'
            with open(filename) as f:
                for line in f:
                    ls = line.strip()
                    if ls:
                        if "#" in ls:
                            continue
                        else:
                            DI_ct_p = []
                            l = re.split('\t',ls)
                            l = [e for e in l if e]
                            #print(l)
                            for s in l:
                                DI_cti = float(s)
                                DI_ct_p.append(DI_cti)
                            DI_ct_p = np.asarray(DI_ct_p)
                        Xe_values_ct.append(DI_ct_p)

            filename = specdist_ct_spectra_lib.path_to_spectra + '_hubble_G_' + str_gamma + '/spectra_'  + case + '_hubble_G_' + str_gamma  + '_Xe_values_ct.txt'
            with open(filename) as f:
                for line in f:
                    ls = line.strip()
                    if ls:
                        if "#" in ls:
                            continue
                        else:
                            DI_ct_p = []
                            l = re.split('\t',ls)
                            l = [e for e in l if e]
                            #print(l)
                            for s in l:
                                DI_cti = float(s)
                                DI_ct_p.append(DI_cti)
                            DI_ct_p = np.asarray(DI_ct_p)
                        Xe_values_no_inj_ct.append(DI_ct_p)
            filename = specdist_ct_spectra_lib.path_to_spectra + '_G_' + str_gamma + '/spectra_' + case + '_G_' + str_gamma  + '_Xe_redshifts_ct.txt'
            with open(filename) as f:
                for line in f:
                    ls = line.strip()
                    if ls:
                        if "#" in ls:
                            continue
                        else:
                            DI_ct_p = []
                            l = re.split('\t',ls)
                            l = [e for e in l if e]
                            #print(l)
                            for s in l:
                                DI_cti = float(s)
                                DI_ct_p.append(DI_cti)
                            DI_ct_p = np.asarray(DI_ct_p)
                        Xe_redshifts_ct.append(DI_ct_p)

            filename = specdist_ct_spectra_lib.path_to_spectra + '_hubble_G_' + str_gamma + '/spectra_'  + case + '_hubble_G_' + str_gamma  + '_Xe_redshifts_ct.txt'
            with open(filename) as f:
                for line in f:
                    ls = line.strip()
                    if ls:
                        if "#" in ls:
                            continue
                        else:
                            DI_ct_p = []
                            l = re.split('\t',ls)
                            l = [e for e in l if e]
                            #print(l)
                            for s in l:
                                DI_cti = float(s)
                                DI_ct_p.append(DI_cti)
                            DI_ct_p = np.asarray(DI_ct_p)
                        Xe_redshifts_no_inj_ct.append(DI_ct_p)

        finj_ct = []
        if case == 'mu_fit':
            filename = specdist_ct_spectra_lib.path_to_spectra + '/spectra_' + specdist_ct_spectra_lib.case_id + '_G_' + str_gamma +'/spectra_spectra_' + specdist_ct_spectra_lib.case_id + '_G_' + str_gamma  + '_finj_ct.txt'
        elif 'xe_history' in case:
            filename = specdist_ct_spectra_lib.path_to_spectra + '_G_' + str_gamma + '/spectra_' + case + '_G_' + str_gamma  + '_finj_ct.txt'
        else:
            filename = specdist_ct_spectra_lib.path_to_spectra + '/spectra_' + specdist_ct_spectra_lib.case_id + '_G_' + str_gamma + '_finj_ct.txt'
        with open(filename) as f:
            for line in f:
                ls = line.strip()
                if ls:
                    if "#" in ls:
                        continue
                    else:
                        finj_ct_p = []
                        l = re.split('\t',ls)
                        l = [e for e in l if e]
                        #print(l)
                        for s in l:
                            finj_cti = float(s)
                            finj_ct_p.append(finj_cti)
                    finj_ct.append(finj_ct_p)

        specdist_ct_spectra_lib.X_2d.append(x_ct)
        specdist_ct_spectra_lib.DI_2d.append(DI_ct)
        specdist_ct_spectra_lib.finj_2d.append(finj_ct)

        if 'xe_history' in case:
            specdist_ct_spectra_lib.Xe_values_2d.append(Xe_values_ct)
            specdist_ct_spectra_lib.Xe_values_no_inj_2d.append(Xe_values_no_inj_ct)

            specdist_ct_spectra_lib.Xe_redshifts_2d.append(Xe_redshifts_ct)
            specdist_ct_spectra_lib.Xe_redshifts_no_inj_2d.append(Xe_redshifts_no_inj_ct)

            #print(np.asarray(Xe_redshifts_ct))
            #print(np.asarray(Xe_values_ct))
            DXe_Xe = []
            DXe_Xe_redshifts = []
            for (p,p_no_inj,z,z_no_inj) in zip(Xe_values_ct,Xe_values_no_inj_ct,Xe_redshifts_ct,Xe_redshifts_no_inj_ct):
                DXe_Xe_p =[]
                DXe_Xe_redshifts_p = []
                #check if nan in any of the arrays:
                Arrays_list = [np.asarray(p),
                              np.asarray(p_no_inj),
                              np.asarray(z),
                              np.asarray(z_no_inj)]
                has_nan = False
                for array in Arrays_list:
                    array_sum = np.sum(array)
                    has_nan += np.isnan(array_sum)
                if np.sum(has_nan):
                    #print('filling with nans')
                    DXe_Xe_p =  np.empty(len(p))
                    DXe_Xe_p[:] = np.nan
                    DXe_Xe_redshifts_p = np.empty(len(p))
                    DXe_Xe_redshifts_p[:] = np.nan
                else:
                    f_Xe = interp1d(np.asarray(z), np.asarray(p))
                    f_Xe_no_inj = interp1d(np.asarray(z_no_inj), np.asarray(p_no_inj))


                    new_z_min  = max(np.min(np.asarray(z_no_inj)),np.min(np.asarray(z)))
                    new_z_max  = min(np.max(np.asarray(z_no_inj)),np.max(np.asarray(z)))

                    Nz = max(len(np.asarray(z)),len(np.asarray(z_no_inj)))
                    new_z = np.logspace(np.log10(new_z_min),np.log10(new_z_max),Nz)

                    new_Xe = f_Xe(new_z)
                    new_Xe_no_inj = f_Xe_no_inj(new_z)
                    DXe_Xe_p = (new_Xe - new_Xe_no_inj)/new_Xe_no_inj
                    DXe_Xe_redshifts_p = new_z
                DXe_Xe.append(DXe_Xe_p)
                DXe_Xe_redshifts.append(DXe_Xe_redshifts_p)

            specdist_ct_spectra_lib.DXe_Xe_2d.append(DXe_Xe)
            specdist_ct_spectra_lib.DXe_Xe_redshifts_2d.append(DXe_Xe_redshifts)






def GetSpectra(Gamma_inj_asked,x_inj_asked,x_asked,specdist_ct_spectra_lib):
    r1 = (Gamma_inj_asked - specdist_ct_spectra_lib.Gamma_inj_min)
    r2 = (specdist_ct_spectra_lib.Gamma_inj_max - Gamma_inj_asked)
    r3 = (x_inj_asked - specdist_ct_spectra_lib.x_inj_min)
    r4 = (specdist_ct_spectra_lib.x_inj_max - x_inj_asked)

    if (r1 < 0) or (r2 < 0) or (r3 < 0) or (r4 < 0):
        #print('filling with nans')
        array_x_asked =  np.empty(len(x_asked))
        array_x_asked[:] = np.nan
        array_S_result = np.empty(len(x_asked))
        array_S_result[:] =  np.nan
        F_gamma_asked_xinj_asked = np.nan

    else:
        X_2d =  specdist_ct_spectra_lib.X_2d
        DI_2d = specdist_ct_spectra_lib.DI_2d
        finj_2d = specdist_ct_spectra_lib.finj_2d

        Gamma_values = specdist_ct_spectra_lib.Gamma_values
        x_inj_values =  specdist_ct_spectra_lib.x_inj_values
        # print(Gamma_values)
        # print(x_inj_values)
        # print(find_nearests(Gamma_values, Gamma_inj_asked))
        id_gamma_low = find_nearests(Gamma_values, Gamma_inj_asked)[0]
        id_gamma_high = find_nearests(Gamma_values, Gamma_inj_asked)[1]

        id_xinj_low = find_nearests(x_inj_values, x_inj_asked)[0]
        id_xinj_high = find_nearests(x_inj_values, x_inj_asked)[1]

        #print(id_gamma_low,id_gamma_high,id_xinj_low,id_xinj_high)

        #DI_2d[id_gamma_low][id_xinj_low]  #DI_2d[id_gamma_low][id_xinj_high]
        #X_2d[id_gamma_low][id_xinj_low]   #X_2d[id_gamma_low][id_xinj_high]

        #DI_2d[id_gamma_high][id_xinj_low]  #DI_2d[id_gamma_high][id_xinj_high]
        #X_2d[id_gamma_high][id_xinj_low]   #X_2d[id_gamma_high][id_xinj_high]
        S = [[X_2d[id_gamma_low][id_xinj_low],DI_2d[id_gamma_low][id_xinj_low]],[X_2d[id_gamma_low][id_xinj_high],DI_2d[id_gamma_low][id_xinj_high]],[X_2d[id_gamma_high][id_xinj_low],DI_2d[id_gamma_high][id_xinj_low]],[X_2d[id_gamma_high][id_xinj_high],DI_2d[id_gamma_high][id_xinj_high]]]
        F = [finj_2d[id_gamma_low][id_xinj_low],finj_2d[id_gamma_low][id_xinj_high],finj_2d[id_gamma_high][id_xinj_low],finj_2d[id_gamma_high][id_xinj_high]]
        dict = {
        "gamma_low": Gamma_values[id_gamma_low],
        "gamma_high": Gamma_values[id_gamma_high],
        "xinj_low": x_inj_values[id_xinj_low],
        "xinj_high": x_inj_values[id_xinj_high],
        "spectra": S,
        "finj": F
        }

        #print(F)

        gamma_low = dict["gamma_low"]
        gamma_high = dict["gamma_high"]
        xinj_low = dict["xinj_low"]
        xinj_high = dict["xinj_high"]
        S = dict["spectra"]
        F = dict["finj"]

        #print(gamma_low)
        #print(gamma_high)

        S_gamma_low_xinj_low = S[0]
        S_gamma_low_xinj_high = S[1]
        S_gamma_high_xinj_low = S[2]
        S_gamma_high_xinj_high = S[3]

        F_gamma_low_xinj_low = F[0][0]
        F_gamma_low_xinj_high = F[1][0]
        F_gamma_high_xinj_low = F[2][0]
        F_gamma_high_xinj_high = F[3][0]

        Gamma_asked = Gamma_inj_asked
        xinj_asked = x_inj_asked

        #check if nan in any of the arrays:
        Arrays_list = [S_gamma_low_xinj_low,
                      S_gamma_low_xinj_high,
                      S_gamma_high_xinj_low,
                      S_gamma_high_xinj_high]
        has_nan = False
        for p in Arrays_list:
            array = p
            array_sum = np.sum(array)
            has_nan += np.isnan(array_sum)
        if np.sum(has_nan):
            #print('filling with nans')
            array_x_asked =  np.empty(len(x_asked))
            array_x_asked[:] = np.nan
            array_S_result = np.empty(len(x_asked))
            array_S_result[:] =  np.nan
            F_gamma_asked_xinj_asked = np.nan
        else:
            nx = int(1e4)


            ############### xinj_low
            new_x_min = np.maximum(np.min(S_gamma_low_xinj_low[0]),np.min(S_gamma_high_xinj_low[0]))
            new_x_max = np.minimum(np.max(S_gamma_low_xinj_low[0]),np.max(S_gamma_high_xinj_low[0]))
            new_x_array = np.logspace(np.log10(new_x_min),np.log10(new_x_max),nx)
            new_x_array = new_x_array[1:-1]


            f_gamma_low = interp1d(S_gamma_low_xinj_low[0], S_gamma_low_xinj_low[1])
            f_gamma_high = interp1d(S_gamma_high_xinj_low[0], S_gamma_high_xinj_low[1])

            new_S_gamma_low = f_gamma_low(new_x_array)
            new_S_gamma_high = f_gamma_high(new_x_array)

            if gamma_low == Gamma_asked:
                w = 1.
            elif gamma_high == Gamma_asked:
                w = 0.
            else:
                #w = (gamma_high - Gamma_asked)/(gamma_high - gamma_low)
                w = (np.log(gamma_high) - np.log(Gamma_asked))/(np.log(gamma_high) - np.log(gamma_low))
            new_S_gamma_asked = w*new_S_gamma_low + (1.-w)*new_S_gamma_high
#             print('xinj _low : w_gamma = %.14e'%w)



            S_gamma_asked_xinj_low = [[],[]]
            S_gamma_asked_xinj_low[0] = new_x_array
            S_gamma_asked_xinj_low[1] = new_S_gamma_asked

            F_gamma_asked_xinj_low = w*F_gamma_low_xinj_low + (1.-w)*F_gamma_high_xinj_low

            ############# xinj_high

            new_x_min = np.maximum(np.min(S_gamma_low_xinj_high[0]),np.min(S_gamma_high_xinj_high[0]))
            new_x_max = np.minimum(np.max(S_gamma_low_xinj_high[0]),np.max(S_gamma_high_xinj_high[0]))
            new_x_array = np.logspace(np.log10(new_x_min),np.log10(new_x_max),nx)
            new_x_array = new_x_array[1:-1]

            f_gamma_low = interp1d(S_gamma_low_xinj_high[0], S_gamma_low_xinj_high[1])
            f_gamma_high = interp1d(S_gamma_high_xinj_high[0], S_gamma_high_xinj_high[1])

            new_S_gamma_low = f_gamma_low(new_x_array)
            new_S_gamma_high = f_gamma_high(new_x_array)

            if gamma_low == Gamma_asked:
                w = 1.
            elif gamma_high == Gamma_asked:
                w = 0.
            else:
                #w = (gamma_high - Gamma_asked)/(gamma_high - gamma_low)
                w = (np.log(gamma_high) - np.log(Gamma_asked))/(np.log(gamma_high) - np.log(gamma_low))
#             #w = (gamma_high - Gamma_asked)/(gamma_high - gamma_low)
#             w = (np.log(gamma_high) - np.log(Gamma_asked))/(np.log(gamma_high) - np.log(gamma_low))
            new_S_gamma_asked = w*new_S_gamma_low + (1.-w)*new_S_gamma_high
#             print('xinj _high : w_gamma = %.14e'%w)

            S_gamma_asked_xinj_high = [[],[]]
            S_gamma_asked_xinj_high[0] = new_x_array
            S_gamma_asked_xinj_high[1] = new_S_gamma_asked

            F_gamma_asked_xinj_high = w*F_gamma_low_xinj_low + (1.-w)*F_gamma_high_xinj_high

            ############# interpolation between xinjs
            new_x_min = np.maximum(np.min(S_gamma_asked_xinj_low[0]),np.min(S_gamma_asked_xinj_high[0]))
            new_x_max = np.minimum(np.max(S_gamma_asked_xinj_low[0]),np.max(S_gamma_asked_xinj_high[0]))
            new_x_array = np.logspace(np.log10(new_x_min),np.log10(new_x_max),nx)
            new_x_array = new_x_array[1:-1]

            f_xinj_low = interp1d(S_gamma_asked_xinj_low[0], S_gamma_asked_xinj_low[1])
            f_xinj_high = interp1d(S_gamma_asked_xinj_high[0], S_gamma_asked_xinj_high[1])

            new_S_xinj_low = f_xinj_low(new_x_array)
            new_S_xinj_high = f_xinj_high(new_x_array)

            #w = (xinj_high - xinj_asked)/(xinj_high - xinj_low)
            w = (np.log(xinj_high) - np.log(xinj_asked))/(np.log(xinj_high) - np.log(xinj_low))
#             print('xinj_high = %.14e'%xinj_high)
#             print('w_xinj = %.14e'%w)

            new_S_xinj_asked = w*new_S_xinj_low + (1.-w)*new_S_xinj_high

            S_gamma_asked_xinj_asked = [[],[]]
            S_gamma_asked_xinj_asked[0] = new_x_array
            S_gamma_asked_xinj_asked[1] = new_S_xinj_asked

            F_gamma_asked_xinj_asked = w*F_gamma_asked_xinj_low + (1.-w)*F_gamma_asked_xinj_high


            f_gamma_asked_xinj_asked = interp1d(S_gamma_asked_xinj_asked[0], S_gamma_asked_xinj_asked[1])
            ########### get spectra at required x values
            bound_x_min = np.min(S_gamma_asked_xinj_asked[0])
            bound_x_max = np.max(S_gamma_asked_xinj_asked[0])

            array_x_asked = np.asarray(x_asked)

            min_x_asked = np.min(array_x_asked)
            max_x_asked = np.max(array_x_asked)

            id_min = 0
            id_max = None
            if min_x_asked < bound_x_min:
                id_min = find_nearests(array_x_asked, bound_x_min)[1]
            if max_x_asked > bound_x_max:
                id_max = find_nearests(array_x_asked, bound_x_max)[0]
            array_x_asked = array_x_asked[id_min:id_max]

            array_S_result = f_gamma_asked_xinj_asked(array_x_asked)

    r_dict = {"x":array_x_asked,
              "DI": array_S_result,
              "finj": F_gamma_asked_xinj_asked}

    return r_dict



def GetXeHistory(Gamma_inj_asked,x_inj_asked,z_asked,specdist_ct_spectra_lib,omega_cdm=0.12,fdm_asked=1e0,get_pca_constraint='yes'):
    x_asked = z_asked
    r1 = (Gamma_inj_asked - specdist_ct_spectra_lib.Gamma_inj_min)
    r2 = (specdist_ct_spectra_lib.Gamma_inj_max - Gamma_inj_asked)
    r3 = (x_inj_asked - specdist_ct_spectra_lib.x_inj_min)
    r4 = (specdist_ct_spectra_lib.x_inj_max - x_inj_asked)
    #print(r1,r2,r3,r4)

    if (r1 < 0) or (r2 < 0) or (r3 < 0) or (r4 < 0):
        #print('filling with nans')
        array_x_asked =  np.empty(len(x_asked))
        array_x_asked[:] = np.nan
        array_S_result = np.empty(len(x_asked))
        array_S_result[:] =  np.nan
        F_gamma_asked_xinj_asked = np.nan

    else:
        X_2d =  specdist_ct_spectra_lib.X_2d
        DI_2d = specdist_ct_spectra_lib.DI_2d
        finj_2d = specdist_ct_spectra_lib.finj_2d

        DXe_Xe_2d = specdist_ct_spectra_lib.DXe_Xe_2d
        Xe_redshifts_2d = specdist_ct_spectra_lib.DXe_Xe_redshifts_2d

        Gamma_values = specdist_ct_spectra_lib.Gamma_values
        x_inj_values =  specdist_ct_spectra_lib.x_inj_values
        # print(Gamma_values)
        # print(x_inj_values)
        # print(find_nearests(Gamma_values, Gamma_inj_asked))
        id_gamma_low = find_nearests(Gamma_values, Gamma_inj_asked)[0]
        id_gamma_high = find_nearests(Gamma_values, Gamma_inj_asked)[1]

        id_xinj_low = find_nearests(x_inj_values, x_inj_asked)[0]
        id_xinj_high = find_nearests(x_inj_values, x_inj_asked)[1]

        #print(id_gamma_low,id_gamma_high,id_xinj_low,id_xinj_high)

        #DI_2d[id_gamma_low][id_xinj_low]  #DI_2d[id_gamma_low][id_xinj_high]
        #X_2d[id_gamma_low][id_xinj_low]   #X_2d[id_gamma_low][id_xinj_high]

        #DI_2d[id_gamma_high][id_xinj_low]  #DI_2d[id_gamma_high][id_xinj_high]
        #X_2d[id_gamma_high][id_xinj_low]   #X_2d[id_gamma_high][id_xinj_high]
        #S = [[X_2d[id_gamma_low][id_xinj_low],DI_2d[id_gamma_low][id_xinj_low]],[X_2d[id_gamma_low][id_xinj_high],DI_2d[id_gamma_low][id_xinj_high]],[X_2d[id_gamma_high][id_xinj_low],DI_2d[id_gamma_high][id_xinj_low]],[X_2d[id_gamma_high][id_xinj_high],DI_2d[id_gamma_high][id_xinj_high]]]
        F = [finj_2d[id_gamma_low][id_xinj_low],finj_2d[id_gamma_low][id_xinj_high],finj_2d[id_gamma_high][id_xinj_low],finj_2d[id_gamma_high][id_xinj_high]]
        DXe_Xe = [[Xe_redshifts_2d[id_gamma_low][id_xinj_low],DXe_Xe_2d[id_gamma_low][id_xinj_low]],
                  [Xe_redshifts_2d[id_gamma_low][id_xinj_high],DXe_Xe_2d[id_gamma_low][id_xinj_high]],
                  [Xe_redshifts_2d[id_gamma_high][id_xinj_low],DXe_Xe_2d[id_gamma_high][id_xinj_low]],
                  [Xe_redshifts_2d[id_gamma_high][id_xinj_high],DXe_Xe_2d[id_gamma_high][id_xinj_high]]]


        dict = {
        "gamma_low": Gamma_values[id_gamma_low],
        "gamma_high": Gamma_values[id_gamma_high],
        "xinj_low": x_inj_values[id_xinj_low],
        "xinj_high": x_inj_values[id_xinj_high],
        "DXe_Xe": DXe_Xe,
        "finj": F
        }

        #print(F)

        gamma_low = dict["gamma_low"]
        gamma_high = dict["gamma_high"]
        xinj_low = dict["xinj_low"]
        xinj_high = dict["xinj_high"]
        S = dict["DXe_Xe"]
        F = dict["finj"]

        #print(gamma_low)
        #print(gamma_high)

        S_gamma_low_xinj_low = S[0]
        S_gamma_low_xinj_high = S[1]
        S_gamma_high_xinj_low = S[2]
        S_gamma_high_xinj_high = S[3]

        F_gamma_low_xinj_low = F[0][0]
        F_gamma_low_xinj_high = F[1][0]
        F_gamma_high_xinj_low = F[2][0]
        F_gamma_high_xinj_high = F[3][0]

        Gamma_asked = Gamma_inj_asked
        xinj_asked = x_inj_asked

        #check if nan in any of the arrays:
        Arrays_list = [S_gamma_low_xinj_low,
                      S_gamma_low_xinj_high,
                      S_gamma_high_xinj_low,
                      S_gamma_high_xinj_high]
        has_nan = False
        for p in Arrays_list:
            array = p
            array_sum = np.sum(array)
            has_nan += np.isnan(array_sum)
        if np.sum(has_nan):
            #print('filling with nans')
            array_x_asked =  np.empty(len(x_asked))
            array_x_asked[:] = np.nan
            array_S_result = np.empty(len(x_asked))
            array_S_result[:] =  np.nan
            F_gamma_asked_xinj_asked = np.nan
        else:
            nx = int(1e4)


            ############### xinj_low
            new_x_min = np.maximum(np.min(S_gamma_low_xinj_low[0]),np.min(S_gamma_high_xinj_low[0]))
            new_x_max = np.minimum(np.max(S_gamma_low_xinj_low[0]),np.max(S_gamma_high_xinj_low[0]))
            new_x_array = np.logspace(np.log10(new_x_min),np.log10(new_x_max),nx)
            new_x_array = new_x_array[1:-1]


            f_gamma_low = interp1d(S_gamma_low_xinj_low[0], S_gamma_low_xinj_low[1])
            f_gamma_high = interp1d(S_gamma_high_xinj_low[0], S_gamma_high_xinj_low[1])

            new_S_gamma_low = f_gamma_low(new_x_array)
            new_S_gamma_high = f_gamma_high(new_x_array)

            if gamma_low == Gamma_asked:
                w = 1.
            elif gamma_high == Gamma_asked:
                w = 0.
            else:
                #w = (gamma_high - Gamma_asked)/(gamma_high - gamma_low)
                w = (np.log(gamma_high) - np.log(Gamma_asked))/(np.log(gamma_high) - np.log(gamma_low))
            new_S_gamma_asked = w*new_S_gamma_low + (1.-w)*new_S_gamma_high
#             print('xinj _low : w_gamma = %.14e'%w)



            S_gamma_asked_xinj_low = [[],[]]
            S_gamma_asked_xinj_low[0] = new_x_array
            S_gamma_asked_xinj_low[1] = new_S_gamma_asked

            F_gamma_asked_xinj_low = w*F_gamma_low_xinj_low + (1.-w)*F_gamma_high_xinj_low

            ############# xinj_high

            new_x_min = np.maximum(np.min(S_gamma_low_xinj_high[0]),np.min(S_gamma_high_xinj_high[0]))
            new_x_max = np.minimum(np.max(S_gamma_low_xinj_high[0]),np.max(S_gamma_high_xinj_high[0]))
            new_x_array = np.logspace(np.log10(new_x_min),np.log10(new_x_max),nx)
            new_x_array = new_x_array[1:-1]

            f_gamma_low = interp1d(S_gamma_low_xinj_high[0], S_gamma_low_xinj_high[1])
            f_gamma_high = interp1d(S_gamma_high_xinj_high[0], S_gamma_high_xinj_high[1])

            new_S_gamma_low = f_gamma_low(new_x_array)
            new_S_gamma_high = f_gamma_high(new_x_array)

            if gamma_low == Gamma_asked:
                w = 1.
            elif gamma_high == Gamma_asked:
                w = 0.
            else:
                #w = (gamma_high - Gamma_asked)/(gamma_high - gamma_low)
                w = (np.log(gamma_high) - np.log(Gamma_asked))/(np.log(gamma_high) - np.log(gamma_low))
#             #w = (gamma_high - Gamma_asked)/(gamma_high - gamma_low)
#             w = (np.log(gamma_high) - np.log(Gamma_asked))/(np.log(gamma_high) - np.log(gamma_low))
            new_S_gamma_asked = w*new_S_gamma_low + (1.-w)*new_S_gamma_high
#             print('xinj _high : w_gamma = %.14e'%w)

            S_gamma_asked_xinj_high = [[],[]]
            S_gamma_asked_xinj_high[0] = new_x_array
            S_gamma_asked_xinj_high[1] = new_S_gamma_asked

            F_gamma_asked_xinj_high = w*F_gamma_low_xinj_low + (1.-w)*F_gamma_high_xinj_high

            ############# interpolation between xinjs
            new_x_min = np.maximum(np.min(S_gamma_asked_xinj_low[0]),np.min(S_gamma_asked_xinj_high[0]))
            new_x_max = np.minimum(np.max(S_gamma_asked_xinj_low[0]),np.max(S_gamma_asked_xinj_high[0]))
            new_x_array = np.logspace(np.log10(new_x_min),np.log10(new_x_max),nx)
            new_x_array = new_x_array[1:-1]

            f_xinj_low = interp1d(S_gamma_asked_xinj_low[0], S_gamma_asked_xinj_low[1])
            f_xinj_high = interp1d(S_gamma_asked_xinj_high[0], S_gamma_asked_xinj_high[1])

            new_S_xinj_low = f_xinj_low(new_x_array)
            new_S_xinj_high = f_xinj_high(new_x_array)

            #w = (xinj_high - xinj_asked)/(xinj_high - xinj_low)
            w = (np.log(xinj_high) - np.log(xinj_asked))/(np.log(xinj_high) - np.log(xinj_low))
#             print('xinj_high = %.14e'%xinj_high)
#             print('w_xinj = %.14e'%w)

            new_S_xinj_asked = w*new_S_xinj_low + (1.-w)*new_S_xinj_high

            S_gamma_asked_xinj_asked = [[],[]]
            S_gamma_asked_xinj_asked[0] = new_x_array
            S_gamma_asked_xinj_asked[1] = new_S_xinj_asked

            F_gamma_asked_xinj_asked = w*F_gamma_asked_xinj_low + (1.-w)*F_gamma_asked_xinj_high


            f_gamma_asked_xinj_asked = interp1d(S_gamma_asked_xinj_asked[0], S_gamma_asked_xinj_asked[1])
            ########### get spectra at required x values
            bound_x_min = np.min(S_gamma_asked_xinj_asked[0])
            bound_x_max = np.max(S_gamma_asked_xinj_asked[0])

            array_x_asked = np.asarray(x_asked)

            min_x_asked = np.min(array_x_asked)
            max_x_asked = np.max(array_x_asked)

            id_min = 0
            id_max = None
            if min_x_asked < bound_x_min:
                id_min = find_nearests(array_x_asked, bound_x_min)[1]
            if max_x_asked > bound_x_max:
                id_max = find_nearests(array_x_asked, bound_x_max)[0]
            array_x_asked = array_x_asked[id_min:id_max]

            array_S_result = f_gamma_asked_xinj_asked(array_x_asked)
    #print(F_gamma_asked_xinj_asked)
    fdm = 1./1.3098e4*F_gamma_asked_xinj_asked*x_inj_asked*(omega_cdm/0.12)**-1
    array_DXe_Xe = array_S_result*fdm_asked/fdm

    #print(fdm)
    if get_pca_constraint == 'yes':
        z1 = specdist_ct_spectra_lib.Xe_PCA_EigenModes['E1']['z']
        E1 = specdist_ct_spectra_lib.Xe_PCA_EigenModes['E1']['values']
        z2 = specdist_ct_spectra_lib.Xe_PCA_EigenModes['E2']['z']
        E2 = specdist_ct_spectra_lib.Xe_PCA_EigenModes['E2']['values']
        z3 = specdist_ct_spectra_lib.Xe_PCA_EigenModes['E3']['z']
        E3 = specdist_ct_spectra_lib.Xe_PCA_EigenModes['E3']['values']

        f_E1 = interp1d(z1, E1)
        f_E2 = interp1d(z2, E2)
        f_E3 = interp1d(z3, E3)
        #print(z)

        f = interp1d(array_x_asked, array_S_result)

        min_z1 = max(np.min(z1),np.min(array_x_asked))
        max_z1 = min(np.max(z1),np.max(array_x_asked))

        min_z2 = max(np.min(z2),np.min(array_x_asked))
        max_z2 = min(np.max(z2),np.max(array_x_asked))

        min_z3 = max(np.min(z3),np.min(array_x_asked))
        max_z3 = min(np.max(z3),np.max(array_x_asked))

        new_z1 = np.linspace(min_z1,max_z1,5000)
        new_z2 = np.linspace(min_z2,max_z2,5000)
        new_z3 = np.linspace(min_z3,max_z3,5000)


        zeta1 = f(new_z1)
        zeta2 = f(new_z2)
        zeta3 = f(new_z3)

        integrand_rho1 = zeta1*f_E1(new_z1)
        integrand_rho2 = zeta2*f_E2(new_z2)
        integrand_rho3 = zeta3*f_E3(new_z3)

        rho1 = np.trapz(integrand_rho1, x=new_z1)/fdm
        rho2 = np.trapz(integrand_rho2, x=new_z2)/fdm
        rho3 = np.trapz(integrand_rho3, x=new_z3)/fdm

        sigma1 = 0.12
        sigma2 = 0.19
        sigma3 = 0.35

        fdm_pca_lim = 2.*np.sqrt(rho1**2/sigma1**2+rho2**2/sigma2**2+rho3**2/sigma3**2)
    else:
        fdm_pca_lim = 0.

    r_dict = {
              "z":array_x_asked,
              "DXe_Xe": array_DXe_Xe,
              "fdm_pca_lim": fdm_pca_lim
              }

    return r_dict
