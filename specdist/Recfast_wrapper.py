from .config import *
from .utils import *



class recfast:
    def __init__(self):
        self.rf_zstart = 2.5e4
        self.rf_zend = 0.0

        self.rf_npts = 10000

        self.rf_T0 = 2.726
        self.rf_N_eff   = 3.046
        self.rf_Yp   = 0.24

        self.rf_Omega_m = 0.26
        self.rf_Omega_b = 0.044

        self.rf_h = 0.71

        self.rf_Recfast_fudge_factor = 0
        self.rf_include_correction_function = 0
        self.rf_A2s1s = 0

        self.rf_Reionization_model = 0
        self.rf_f_ann = 0.0e-24

        self.rf_f_dec = 1e-5
        self.rf_Gamma_dec = 1.e-15
        self.rf_xinj0 = 1.e-4
        self.rf_decay_scenario = 'heating'

        self.rf_B0 = 0.0
        self.rf_nB = -2.9

        self.rf_include_turbulent_decay = 1
        self.rf_include_ambipolar_diffusion = 0
        self.rf_Paoletti_Finelli_Lorentz    = 0

        self.rf_add_collisions = 0

        self.rf_alp_alp_ref       = 1.0
        self.rf_me_me_ref       = 1.0
        self.rf_power_for_1pz_a = 0.0
        self.rf_variation_mode = 0

        self.rf_verbosity_level_Recfast = 1
        self.rf_path_for_output =  './outputs/'
        self.rf_addition_to_name_for_output = '.tmp.dat'

        self.path_to_rf_param_file = path_to_recfast + '/runfiles/'
        self.tmp_dir_name = 'tmp'
        self.save_dir_name = 'tmp'

        #######################


    def create_tmp_dir_to_store_full_rf_outputs(self):
        self.tmp_dir_name = self.save_dir_name
        self.path_to_rf_tmp_dir = self.path_to_rf_param_file + self.tmp_dir_name
        subprocess.call(['rm','-rf',self.path_to_rf_param_file+self.tmp_dir_name])
        subprocess.call(['mkdir',self.path_to_rf_param_file+self.tmp_dir_name])





    def compute_recfast(self,index_pval=0,**params_values_dict):

        p_dict = params_values_dict
        subprocess.call(['mkdir',self.path_to_rf_tmp_dir+'/tmp_'+str(index_pval)])
        p_dict['path for output'] = self.path_to_rf_tmp_dir+'/tmp_'+str(index_pval) + '/'
        with open(self.path_to_rf_tmp_dir+'/tmp_'+str(index_pval)+'/tmp.ini', 'w') as f:
            for k, v in p_dict.items():
                f.write(str(k) + ' = '+ str(v) + '\n')
        f.close()
        subprocess.call([path_to_recfast+'/Recfast++',self.path_to_rf_tmp_dir+'/tmp_'+str(index_pval)+'/tmp.ini'])
        r_dict = {}

        R = np.loadtxt(p_dict['path for output']+'Xe_Recfast++'+self.root_name+'tmp.dat')
        try:
            r_dict['z'] = R[:,0]
            r_dict['Xe'] = R[:,1]
        except IndexError:
            a = np.empty(1)
            a[:] = np.nan
            r_dict['z'] = a
            r_dict['Xe'] = a
        #print(r_dict)
        return r_dict

    def compute_recfast_parallel(self,index_pval,param_values_array,param_name,**kwargs):

        # dict_for_fisher = kwargs.get('dict_for_fisher')
        # sd_lib_for_fisher = kwargs.get('sd_lib_for_fisher', None)

        p_val = param_values_array[index_pval]
        params_values_dict = self.load_parameter_file()
        params_values_dict[param_name] = p_val
        # if float(params_values_dict['pi_f_dm']) != 0 and float(params_values_dict['photon injection f_dec']) == 0:
        #     params_values_dict['photon injection f_dec'] = 1.3098e4*float(params_values_dict['pi_f_dm'])/float(params_values_dict['photon injection x_dec'])*(self.ct_omega_cdm/0.12)*(float(params_values_dict['T0'])/2.726)**-4
        # if params_values_dict['pi_finj_from_fisher'] == 'yes':
        #     f_dm_fisher = pi_run_fisher_constraints([float(params_values_dict['photon injection Gamma_dec'])],[float(params_values_dict['photon injection x_dec'])],sd_lib_for_fisher,**dict_for_fisher)
        #     params_values_dict['photon injection f_dec'] = f_dm_fisher['curves'][0]['finj'][0]
        #     print('finj_fisher = %e'%params_values_dict['photon injection f_dec'])
        r_dict = self.compute_recfast(index_pval=index_pval,**params_values_dict)
        dict_rf_results = r_dict
        dict_param_values = {}
        dict_param_values[param_name] = p_val
        r_dict = {**dict_param_values,**dict_rf_results}
        return r_dict

    def run_recfast_parallel(self,**args):
        self.create_tmp_dir_to_store_full_rf_outputs()
        startTime = datetime.now()
        pool = multiprocessing.Pool()
        if type(args['param_values_array'])== float or type(args['param_values_array'])== int:
            array_args = [args['param_values_array']]
        else:
            array_args = args['param_values_array']

        fn=functools.partial(self.compute_recfast_parallel,param_values_array=array_args,param_name=args['param_name'])
        #print(len(*param_values_array))
        results = pool.map(fn,range(np.size(np.asarray(args['param_values_array']))))
        pool.close()
        #self.clear()
        #if self.ct_pi_redshift_evolution_mode==0:
        try:
            if args['save_recfast_results']=='yes':
                subprocess.call(['rm','-rf',path_to_recfast_results+'/'+self.save_dir_name])
                subprocess.call(['mkdir',path_to_recfast_results+'/'+self.save_dir_name])
                z_rf = []
                Xe_rf = []
                # if self.ct_include_pi == 1:
                #     finj_ct = []
                # if self.save_Xe == 'yes' and self.ct_evolve_Xe != 0 :
                #     Xe_values_ct = []
                #     Xe_redshifts_ct = []
                for ip in range(len(results)):
                    z_rf.append(results[ip]['z'])
                    Xe_rf.append(results[ip]['Xe'])
                    # if self.ct_include_pi == 1:
                    #     finj_ct.append(results[ip]['finj'])
                    # if self.save_Xe == 'yes' and self.ct_evolve_Xe != 0 :
                    #     Xe_values_ct.append(results[ip]['Xe_values'])
                    #     Xe_redshifts_ct.append(results[ip]['Xe_redshifts'])

                str_param = 'p'
                if args['param_name'] == 'photon injection x_dec':
                    str_param = 'xinj'
                with open(path_to_recfast_results+'/'+self.save_dir_name + '/' + self.save_dir_name  + '_'+str_param+'_recfast.txt', 'w') as f:
                    f.write("# arrays of %s for recfast\n"%args['param_name'])
                    for row in array_args:
                        np.savetxt(f,[row],fmt="%.3e",delimiter='\t')
                f.close()
                # if self.ct_include_pi == 1:
                #     with open(path_to_ct_spectra_results+'/'+self.save_dir_name + '/' + self.save_dir_name  + '_finj_recfast.txt', 'w') as f:
                #         f.write("# arrays of finj values for CT spectra\n")
                #         for row in finj_ct:
                #             np.savetxt(f,[row],fmt="%.3e",delimiter='\t')
                #     f.close()
                with open(path_to_recfast_results+'/'+self.save_dir_name + '/' + self.save_dir_name  + '_Xe_recfast.txt', 'w') as f:
                    f.write("# arrays of Xe values\n")
                    for row in Xe_rf:
                        np.savetxt(f,[row],fmt="%.8e",delimiter='\t')
                f.close()
                with open(path_to_recfast_results+'/'+self.save_dir_name + '/' + self.save_dir_name  + '_z_recfast.txt', 'w') as f:
                    f.write("# arrays of redshift values\n")
                    for row in z_rf:
                        np.savetxt(f,[row],fmt="%.3e",delimiter='\t')
                f.close()
                # if self.save_Xe == 'yes' and self.ct_evolve_Xe != 0 :
                #     with open(path_to_ct_spectra_results+'/'+self.save_dir_name + '/spectra_' + self.save_dir_name  + '_Xe_redshifts_ct.txt', 'w') as f:
                #         f.write("# arrays of redshift values for free electron fraction Xe\n")
                #         for row in Xe_redshifts_ct:
                #             np.savetxt(f,[row],fmt="%.3e",delimiter='\t')
                #     f.close()
                #     with open(path_to_ct_spectra_results+'/'+self.save_dir_name + '/spectra_' + self.save_dir_name  + '_Xe_values_ct.txt', 'w') as f:
                #         f.write("# arrays of Xe values for free electron fraction Xe\n")
                #         for row in Xe_values_ct:
                #             np.savetxt(f,[row],fmt="%.8e",delimiter='\t')
                #     f.close()
        except KeyError:
            pass



        return results


    def load_parameter_file(self):
        #load template parameter file into dictionnary
        p_dict = {}
        with open(self.path_to_rf_param_file+"recfast-parameters.ini") as f:
            for line in f:
                x = line.strip()
                if x:
                    if not x.startswith("#"):
                        l = re.split(r'[=#]',x)
                        (key, val) = (l[0].strip(),l[1].strip())
                        p_dict[key] = val
        f.close()

        p_dict['zstart'] = self.rf_zstart
        p_dict['zend'] = self.rf_zend
        p_dict['npts'] = self.rf_npts

        p_dict['T0'] = self.rf_T0
        p_dict['Yp'] = self.rf_Yp
        p_dict['N_eff'] = self.rf_N_eff

        p_dict['Omega_m'] = self.rf_Omega_m
        p_dict['Omega_b'] = self.rf_Omega_b

        p_dict['h100'] = self.rf_h

        p_dict['f_dec'] = self.rf_f_dec
        p_dict['Gamma_dec'] = self.rf_Gamma_dec
        p_dict['xinj0'] = self.rf_xinj0
        p_dict['decay scenario'] = self.rf_decay_scenario
        p_dict['include correction function'] = self.rf_include_correction_function

        p_dict['Reionization model'] = self.rf_Reionization_model
        p_dict['add collisions'] = self.rf_add_collisions

        p_dict['path for output'] = '/outputs/'
        p_dict['addition to name for output'] = '.tmp.dat'
        p_dict['verbosity level Recfast++'] = self.rf_verbosity_level_Recfast
        if self.rf_include_correction_function == 0 and self.rf_Reionization_model == 0:
            self.root_name = '.DM_decay_heating.Tr.'
        elif self.rf_include_correction_function == 1 and self.rf_Reionization_model == 0:
            self.root_name = '.DM_decay_heating.Rec_corrs_CT2010.Tr.'
        elif self.rf_include_correction_function == 0 and self.rf_Reionization_model == 1:
            self.root_name = '.DM_decay_heating.Tr.reion_VP.'
        elif self.rf_include_correction_function == 1 and self.rf_Reionization_model == 1:
            self.root_name = '.DM_decay_heating.Rec_corrs_CT2010.Tr.reion_VP.'
        ##############
        return p_dict

    def clear(self):
        subprocess.call(['rm','-rf',self.path_to_rf_tmp_dir])
