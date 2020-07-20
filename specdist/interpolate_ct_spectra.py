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



        finj_ct = []
        if case == 'mu_fit':
            filename = specdist_ct_spectra_lib.path_to_spectra + '/spectra_' + specdist_ct_spectra_lib.case_id + '_G_' + str_gamma +'/spectra_spectra_' + specdist_ct_spectra_lib.case_id + '_G_' + str_gamma  + '_finj_ct.txt'
        elif 'xe_history' in case:
            specdist_ct_spectra_lib.path_to_spectra + '_G_' + str_gamma + '/spectra_' + case + '_G_' + str_gamma  + '_finj_ct.txt'
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



def GetSpectra(Gamma_inj_asked,x_inj_asked,x_asked,specdist_ct_spectra_lib):
    r1 = (Gamma_inj_asked - specdist_ct_spectra_lib.Gamma_inj_min)
    r2 = (specdist_ct_spectra_lib.Gamma_inj_max - Gamma_inj_asked)
    r3 = (x_inj_asked - specdist_ct_spectra_lib.x_inj_min)
    r4 = (specdist_ct_spectra_lib.x_inj_max - x_inj_asked)

    if (r1 < 0) or (r2 < 0) or (r3 < 0) or (r4 < 0):
        print('filling with nans')
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
        if has_nan:
            print('filling with nans')
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

            #w = (gamma_high - Gamma_asked)/(gamma_high - gamma_low)
            w = (np.log(gamma_high) - np.log(Gamma_asked))/(np.log(gamma_high) - np.log(gamma_low))
            new_S_gamma_asked = w*new_S_gamma_low + (1.-w)*new_S_gamma_high



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

            #w = (gamma_high - Gamma_asked)/(gamma_high - gamma_low)
            w = (np.log(gamma_high) - np.log(Gamma_asked))/(np.log(gamma_high) - np.log(gamma_low))
            new_S_gamma_asked = w*new_S_gamma_low + (1.-w)*new_S_gamma_high

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
