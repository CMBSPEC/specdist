import subprocess

# This path needs to be adjsuted (TBD: set this automatically)
# path_to_sd_projects = "/Users/boris/Work/SPECTRAL-DISTORTIONS/"
path_to_sd_projects = "/scratch/nas_chluba/specdist/bolliet/specdist_ml/"

#path to the cosmotherm database
path_to_ct_database = path_to_sd_projects + "specdist/specdist/data/ct_database/"


#the path to the cosmotherm binary file
path_to_cosmotherm = path_to_sd_projects + "cosmotherm.rel_corr"

#the path to save the output from cosmotherm
path_to_ct_spectra_results =  path_to_sd_projects + "specdist/specdist/ct_spectra"

subprocess.call(['mkdir','-p',path_to_ct_spectra_results])


#the path to the Recfast binary file
path_to_recfast = path_to_sd_projects + "recfast-.vX"

#the path to save the output from cosmotherm
path_to_recfast_results =  path_to_sd_projects + "specdist/specdist/recfast_outputs"

subprocess.call(['mkdir','-p',path_to_recfast_results])
