from .config import *
from .utils import *


class cosmotherm:
    ct_x_dec = 1e-6
    ct_Drho_rho_dec = 3e-5
    ct_npts = 3000
    ct_zstart = 5.e6
    ct_zend = 1e-2
    ct_Gamma_dec = 1e-9
    ct_verbose = 2
    ct_lyc = 0
    ct_pi_redshift_evolution_mode = 0
    ct_include_pi = 1 #photon injection case

    ct_Omega_m = 0.26
    ct_Omega_b = 0.044
    ct_T0 = 2.726
    ct_h = 0.71
    ct_N_eff   = 3.046
    ct_Yp   = 0.24

    path_to_ct_param_file = path_to_cosmotherm + '/runfiles/'

    def __init__(self):
        subprocess.call(['rm','-rf',self.path_to_ct_param_file+'tmp'])
        subprocess.call(['mkdir',self.path_to_ct_param_file+'tmp'])



    def compute_specdist(self,index_pval=0,**params_values_dict):
        p_dict = params_values_dict
        subprocess.call(['rm','-rf',self.path_to_ct_param_file+'tmp'+str(index_pval)])
        subprocess.call(['mkdir',self.path_to_ct_param_file+'tmp/tmp_'+str(index_pval)])
        p_dict['path for output'] = self.path_to_ct_param_file+'tmp/tmp_'+str(index_pval) + '/'
        with open(self.path_to_ct_param_file+'tmp/tmp_'+str(index_pval)+'/tmp.ini', 'w') as f:
            for k, v in p_dict.items():
                f.write(str(k) + ' = '+ str(v) + '\n')
        f.close()
        subprocess.call([path_to_cosmotherm+'/CosmoTherm',self.path_to_ct_param_file+'tmp/tmp_'+str(index_pval)+'/tmp.ini'])
        r_dict = {}
        if (self.ct_pi_redshift_evolution_mode == 1):
            R = np.loadtxt(p_dict['path for output']+'pi_finj_calc.txt')
            r_dict["z"] = R[:,0]
            r_dict["dDrho_rhodt_rel"] = R[:,2]
            with open(p_dict['path for output']+'pi_finj_calc.txt') as f:
                line = f.readline()
                x = line.strip()
                l = re.split(r'[=#\t]',x)
                l[:] = [e.strip() for e in l if e]
                l_dict_keys = l[0::2]
                l_dict_values = l[1::2]
                l_dict = {}
                for k,v in zip(l_dict_keys,l_dict_values):
                    l_dict[k]=float(v)
            r_dict = {**r_dict, **l_dict}
        else:
            R = np.loadtxt(p_dict['path for output']+'Dn.cooling'+self.root_name+'PDE_ODE.tmp.dat')
            r_dict['x'] = R[:,0]
            r_dict['DI'] = R[:,5]
        return r_dict

    def compute_specdist_parallel(self,index_pval,param_values_array,param_name):
        p_val = param_values_array[index_pval]
        params_values_dict = self.load_parameter_file()
        params_values_dict[param_name] = p_val
        r_dict = self.compute_specdist(index_pval=index_pval,**params_values_dict)
        dict_ct_results = r_dict
        dict_param_values = {}
        dict_param_values[param_name] = p_val
        r_dict = {**dict_param_values,**dict_ct_results}
        return r_dict

    def run_cosmotherm_parallel(self,**args):
        startTime = datetime.now()
        pool = multiprocessing.Pool()
        fn=functools.partial(self.compute_specdist_parallel,param_values_array=args['param_values_array'],param_name=args['param_name'])
        #print(len(*param_values_array))
        results = pool.map(fn,range(len(args['param_values_array'])))
        pool.close()
        return results


    def load_parameter_file(self):
        #load template parameter file into dictionnary
        p_dict = {}
        with open(self.path_to_ct_param_file+"parameters.ini") as f:
            for line in f:
                x = line.strip()
                if x:
                    if not x.startswith("#"):
                        l = re.split(r'[=#]',x)
                        (key, val) = (l[0].strip(),l[1].strip())
                        p_dict[key] = val
        f.close()

        p_dict['path for output'] = '/outputs/'
        p_dict['addition to filename at end'] = '.tmp.dat'
        p_dict['N_eff'] = self.ct_N_eff
        p_dict['Yp'] = self.ct_Yp
        p_dict['include photon injection from decaying particle'] = self.ct_include_pi
        if (p_dict['include photon injection from decaying particle'] == 1):
            self.root_name = '.photon_inj.'
        else:
            self.root_name = '.'
        # <!> mX_dec_in_eV must not be passed for x_dec to be read
        p_dict['photon injection x_dec'] = self.ct_x_dec
        p_dict['photon injection Drho_rho_dec'] = self.ct_Drho_rho_dec
        p_dict['npts'] = self.ct_npts
        p_dict['zstart'] = self.ct_zstart
        p_dict['zend'] = self.ct_zend
        p_dict['verbosity level CosmoTherm'] = self.ct_verbose
        p_dict['photon injection Gamma_dec'] = self.ct_Gamma_dec
        p_dict['include Lyc absorption'] = self.ct_lyc
        p_dict['photon injection redshift evolution'] = self.ct_pi_redshift_evolution_mode
        return p_dict

    def clear(self):
        subprocess.call(['rm','-rf',self.path_to_ct_param_file+'tmp'])
