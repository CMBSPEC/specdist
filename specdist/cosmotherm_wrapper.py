from .config import *
from .utils import *
from .pispec_run_fisher_constraints import *


class cosmotherm:
    def __init__(self):
        self.ct_x_dec = 1e-6
        self.ct_Drho_rho_dec = 3e-5
        self.ct_photon_injection_f_dec = 0.
        self.ct_npts = 3000
        self.ct_zstart = 5.e6
        self.ct_zend = 1e-2
        self.ct_zlate = -1
        self.ct_Gamma_dec = 1e-9
        self.ct_verbose = 2
        self.ct_lyc = 0
        self.ct_evolve_Xe = 0
        self.ct_pi_redshift_evolution_mode = 0
        self.ct_pi_finj_mode = 0
        self.ct_pi_stim = 0
        self.ct_include_collisions = 0
        self.ct_include_pi = 1 #photon injection case
        self.ct_reionisation_model = 0
        self.ct_emission_absorption_mode = 0
        self.ct_Omega_m = 0.26
        self.ct_Omega_b = 0.044
        self.ct_h = 0.71
        self.ct_omega_cdm = (self.ct_Omega_m - self.ct_Omega_b)*self.ct_h**2.
        self.ct_T0 = 2.726
        self.ct_N_eff   = 3.046
        self.ct_Yp   = 0.24
        self.save_dir_name = 'tmp'
        self.save_Xe = 'no'
        self.ct_fdm = 0
        self.get_finj = 0
        self.ct_pi_energy_norm = 0
        self.ct_pi_finj_from_fisher = 'no'

        self.path_to_ct_param_file = path_to_cosmotherm + '/runfiles/'
        self.tmp_dir_name = 'tmp'



    def create_tmp_dir_to_store_full_ct_outputs(self):
        self.tmp_dir_name = self.save_dir_name
        self.path_to_ct_tmp_dir = self.path_to_ct_param_file + self.tmp_dir_name
        subprocess.call(['rm','-rf',self.path_to_ct_param_file+self.tmp_dir_name])
        subprocess.call(['mkdir',self.path_to_ct_param_file+self.tmp_dir_name])





    def compute_specdist(self,index_pval=0,**params_values_dict):

        p_dict = params_values_dict
        subprocess.call(['mkdir',self.path_to_ct_tmp_dir+'/tmp_'+str(index_pval)])
        p_dict['path for output'] = self.path_to_ct_tmp_dir+'/tmp_'+str(index_pval) + '/'
        with open(self.path_to_ct_tmp_dir+'/tmp_'+str(index_pval)+'/tmp.ini', 'w') as f:
            for k, v in p_dict.items():
                f.write(str(k) + ' = '+ str(v) + '\n')
        f.close()
        subprocess.call([path_to_cosmotherm+'/CosmoTherm',self.path_to_ct_tmp_dir+'/tmp_'+str(index_pval)+'/tmp.ini'])
        r_dict = {}
        if(self.get_finj==1):
            R = np.loadtxt(p_dict['path for output']+'pi_finj_boris_jens.txt')
            r_dict['finj_boris'] = R[0]
            r_dict['finj_jens'] = R[1]
        else:
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
                try:
                    r_dict['x'] = R[:,0]
                    r_dict['DI'] = R[:,5]
                    if self.save_Xe == 'yes' and self.ct_evolve_Xe != 0 :
                        R = np.loadtxt(p_dict['path for output']+'Xe_Xp_etc.cooling'+self.root_name+'PDE_ODE.tmp.dat')
                        r_dict['Xe_redshifts'] = R[:,0]
                        r_dict['Xe_values'] = R[:,6]
                    if self.ct_include_pi == 1:
                        f = open(p_dict['path for output']+'/parameter_info.cooling'+self.root_name+'PDE_ODE.tmp.dat')
                        lines = f.readlines()
                        for line in lines:
                            if 'finj =' in line:
                                for t in line.split():
                                    try:
                                        finj = float(line.split()[2])
                                    except ValueError:
                                        print('error for process %d, f_inj not found'%index_pval)
                                        pass
                        f.close()
                        r_dict['finj'] = finj
                except IndexError:
                    a = np.empty(1)
                    a[:] = np.nan
                    r_dict['x'] = a
                    r_dict['DI'] = a
                    if self.save_Xe == 'yes' and self.ct_evolve_Xe != 0 :
                        r_dict['Xe_redshifts'] = a
                        r_dict['Xe_values'] = a
                    if self.ct_include_pi == 1:
                        r_dict['finj'] = a[0]

        #print(r_dict)
        return r_dict

    def compute_specdist_parallel(self,index_pval,param_values_array,param_name,**kwargs):

        dict_for_fisher = kwargs.get('dict_for_fisher')
        sd_lib_for_fisher = kwargs.get('sd_lib_for_fisher', None)

        p_val = param_values_array[index_pval]
        params_values_dict = self.load_parameter_file()
        params_values_dict[param_name] = p_val
        if float(params_values_dict['pi_f_dm']) != 0 and float(params_values_dict['photon injection f_dec']) == 0:
            params_values_dict['photon injection f_dec'] = 1.3098e4*float(params_values_dict['pi_f_dm'])/float(params_values_dict['photon injection x_dec'])*(self.ct_omega_cdm/0.12)*(float(params_values_dict['T0'])/2.726)**-4
        if params_values_dict['pi_finj_from_fisher'] == 'yes':
            f_dm_fisher = pi_run_fisher_constraints([float(params_values_dict['photon injection Gamma_dec'])],[float(params_values_dict['photon injection x_dec'])],sd_lib_for_fisher,**dict_for_fisher)
            params_values_dict['photon injection f_dec'] = f_dm_fisher['curves'][0]['finj'][0]
            print('finj_fisher = %e'%params_values_dict['photon injection f_dec'])
        r_dict = self.compute_specdist(index_pval=index_pval,**params_values_dict)
        dict_ct_results = r_dict
        dict_param_values = {}
        dict_param_values[param_name] = p_val
        r_dict = {**dict_param_values,**dict_ct_results}
        return r_dict

    def run_cosmotherm_parallel(self,**args):
        dict_for_fisher = args.get('dict_for_fisher')
        sd_lib_for_fisher = args.get('sd_lib_for_fisher', None)
        args['dict_for_fisher'] = dict_for_fisher
        args['sd_lib_for_fisher'] = sd_lib_for_fisher

        self.create_tmp_dir_to_store_full_ct_outputs()
        startTime = datetime.now()
        pool = multiprocessing.Pool()
        if type(args['param_values_array'])== float or type(args['param_values_array'])== int:
            array_args = [args['param_values_array']]
        else:
            array_args = args['param_values_array']

        fn=functools.partial(self.compute_specdist_parallel,param_values_array=array_args,param_name=args['param_name'],dict_for_fisher=args['dict_for_fisher'],sd_lib_for_fisher=args['sd_lib_for_fisher'])
        #print(len(*param_values_array))
        results = pool.map(fn,range(np.size(np.asarray(args['param_values_array']))))
        pool.close()
        #self.clear()
        if self.ct_pi_redshift_evolution_mode==0:
            try:
                if args['save_spectra']=='yes':
                    subprocess.call(['rm','-rf',path_to_ct_spectra_results+'/'+self.save_dir_name])
                    subprocess.call(['mkdir',path_to_ct_spectra_results+'/'+self.save_dir_name])
                    x_ct = []
                    DI_ct = []
                    if self.ct_include_pi == 1:
                        finj_ct = []
                    if self.save_Xe == 'yes' and self.ct_evolve_Xe != 0 :
                        Xe_values_ct = []
                        Xe_redshifts_ct = []
                    for ip in range(len(results)):
                        x_ct.append(results[ip]['x'])
                        DI_ct.append(results[ip]['DI'])
                        if self.ct_include_pi == 1:
                            finj_ct.append(results[ip]['finj'])
                        if self.save_Xe == 'yes' and self.ct_evolve_Xe != 0 :
                            Xe_values_ct.append(results[ip]['Xe_values'])
                            Xe_redshifts_ct.append(results[ip]['Xe_redshifts'])

                    str_param = 'p'
                    if args['param_name'] == 'photon injection x_dec':
                        str_param = 'xinj'
                    with open(path_to_ct_spectra_results+'/'+self.save_dir_name + '/spectra_' + self.save_dir_name  + '_'+str_param+'_ct.txt', 'w') as f:
                        f.write("# arrays of %s for CT spectra\n"%args['param_name'])
                        for row in array_args:
                            np.savetxt(f,[row],fmt="%.3e",delimiter='\t')
                    f.close()
                    if self.ct_include_pi == 1:
                        with open(path_to_ct_spectra_results+'/'+self.save_dir_name + '/spectra_' + self.save_dir_name  + '_finj_ct.txt', 'w') as f:
                            f.write("# arrays of finj values for CT spectra\n")
                            for row in finj_ct:
                                np.savetxt(f,[row],fmt="%.3e",delimiter='\t')
                        f.close()
                    with open(path_to_ct_spectra_results+'/'+self.save_dir_name + '/spectra_' + self.save_dir_name  + '_DI_ct.txt', 'w') as f:
                        f.write("# arrays of DI values for CT spectra\n")
                        for row in DI_ct:
                            np.savetxt(f,[row],fmt="%.8e",delimiter='\t')
                    f.close()
                    with open(path_to_ct_spectra_results+'/'+self.save_dir_name + '/spectra_' + self.save_dir_name  + '_x_ct.txt', 'w') as f:
                        f.write("# arrays of x=hnu/kT values for CT spectra\n")
                        for row in x_ct:
                            np.savetxt(f,[row],fmt="%.3e",delimiter='\t')
                    f.close()
                    if self.save_Xe == 'yes' and self.ct_evolve_Xe != 0 :
                        with open(path_to_ct_spectra_results+'/'+self.save_dir_name + '/spectra_' + self.save_dir_name  + '_Xe_redshifts_ct.txt', 'w') as f:
                            f.write("# arrays of redshift values for free electron fraction Xe\n")
                            for row in Xe_redshifts_ct:
                                np.savetxt(f,[row],fmt="%.3e",delimiter='\t')
                        f.close()
                        with open(path_to_ct_spectra_results+'/'+self.save_dir_name + '/spectra_' + self.save_dir_name  + '_Xe_values_ct.txt', 'w') as f:
                            f.write("# arrays of Xe values for free electron fraction Xe\n")
                            for row in Xe_values_ct:
                                np.savetxt(f,[row],fmt="%.8e",delimiter='\t')
                        f.close()
            except KeyError:
                pass



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
        p_dict['Omega_m'] = self.ct_Omega_m
        p_dict['Omega_b'] = self.ct_Omega_b
        p_dict['h100'] = self.ct_h
        p_dict['T0'] = self.ct_T0
        p_dict['include photon injection from decaying particle'] = self.ct_include_pi
        if (p_dict['include photon injection from decaying particle'] == 1):
            self.root_name = '.photon_inj.'
        else:
            self.root_name = '.'
        # <!> mX_dec_in_eV must not be passed for x_dec to be read
        p_dict['photon injection x_dec'] = self.ct_x_dec
        p_dict['photon injection Drho_rho_dec'] = self.ct_Drho_rho_dec
        # if self.ct_fdm != 0 and self.ct_photon_injection_f_dec == 0:
        #     self.ct_photon_injection_f_dec = 1.3098e4*self.ct_f_dm/self.ct_x_dec*(self.ct_omega_cdm/0.12)*(self.ct_T0/2.726)**-4
        p_dict['photon injection f_dec'] = self.ct_photon_injection_f_dec
        p_dict['npts'] = self.ct_npts
        p_dict['zstart'] = self.ct_zstart
        p_dict['zend'] = self.ct_zend
        p_dict['zlate'] = self.ct_zlate
        p_dict['verbosity level CosmoTherm'] = self.ct_verbose
        p_dict['photon injection Gamma_dec'] = self.ct_Gamma_dec
        p_dict['include Lyc absorption'] = self.ct_lyc
        p_dict['evolve Xe'] = self.ct_evolve_Xe
        p_dict['Reionization model'] = self.ct_reionisation_model
        p_dict['photon injection redshift evolution'] = self.ct_pi_redshift_evolution_mode
        p_dict['photon injection finj mode'] = self.ct_pi_finj_mode
        p_dict['photon injection energy norm'] = self.ct_pi_energy_norm
        p_dict['photon injection stimulated'] = self.ct_pi_stim
        p_dict['emission/absorption mode'] = self.ct_emission_absorption_mode
        p_dict['pi_f_dm'] = self.ct_fdm
        p_dict['pi_finj_from_fisher'] = self.ct_pi_finj_from_fisher
        p_dict['include collisions'] = self.ct_include_collisions
        return p_dict

    def compute_no_injection_spectrum_and_Xe_history(self,**kwargs):
        r_dict = {}
        self.ct_Drho_rho_dec = 1e-300

        Nx_no_inj = kwargs.get('Nx_no_inj', 5e3)
        Nz_no_inj = kwargs.get('Nz_no_inj', 5e3)
        xmin = kwargs.get('xmin_no_inj', 1e-10)
        xmax = kwargs.get('xmax_no_inj', 1e10)

        R = self.run_cosmotherm_parallel(**kwargs)

        f_DI_no_inj = interp1d(R[0]['x'],R[0]['DI'],bounds_error=False,fill_value=1e-300)

        new_x = np.logspace(np.log10(xmin),np.log10(xmax),Nx_no_inj)
        new_DI_no_inj = f_DI_no_inj(new_x)

        r_dict['x'] = new_x
        r_dict['DI'] = new_DI_no_inj

        if self.save_Xe == 'yes' and self.ct_evolve_Xe != 0 :
            f_Xe_no_inj = interp1d(R[0]['Xe_redshifts'],R[0]['Xe_values'],bounds_error=False,fill_value=1e-300)
            new_Xe_redshifts_ct = np.logspace(np.log10(self.ct_zend),np.log10(self.ct_zstart),Nz_no_inj)
            new_Xe_values_ct = f_Xe_no_inj(new_Xe_redshifts_ct)
            r_dict['Xe_redshifts'] = new_Xe_redshifts_ct
            r_dict['Xe_values'] = new_Xe_values_ct

        if kwargs['save_spectra']=='yes':
            with open(path_to_ct_spectra_results+'/'+self.save_dir_name + '/spectra_' + self.save_dir_name  + '_fine_grid_DI_ct.txt', 'w') as f:
                f.write("# arrays of DI values for CT spectra\n")
                for row in [new_DI_no_inj]:
                    np.savetxt(f,[row],fmt="%.8e",delimiter='\t')
            f.close()
            with open(path_to_ct_spectra_results+'/'+self.save_dir_name + '/spectra_' + self.save_dir_name  + '_fine_grid_x_ct.txt', 'w') as f:
                f.write("# arrays of x values for CT spectra\n")
                for row in [new_x]:
                    np.savetxt(f,[row],fmt="%.8e",delimiter='\t')
            f.close()
            if self.save_Xe == 'yes' and self.ct_evolve_Xe != 0 :
                with open(path_to_ct_spectra_results+'/'+self.save_dir_name + '/spectra_' + self.save_dir_name  + '_fine_grid_Xe_redshifts_ct.txt', 'w') as f:
                    f.write("# arrays of redshift values for free electron fraction Xe\n")
                    for row in [new_Xe_redshifts_ct]:
                        np.savetxt(f,[row],fmt="%.3e",delimiter='\t')
                f.close()
                with open(path_to_ct_spectra_results+'/'+self.save_dir_name + '/spectra_' + self.save_dir_name  + '_fine_grid_Xe_values_ct.txt', 'w') as f:
                    f.write("# arrays of Xe values for free electron fraction Xe\n")
                    for row in [new_Xe_values_ct]:
                        np.savetxt(f,[row],fmt="%.8e",delimiter='\t')
                f.close()
        return r_dict




    def clear(self):
        subprocess.call(['rm','-rf',self.path_to_ct_tmp_dir])
