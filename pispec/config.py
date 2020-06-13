import subprocess

#the path to the cosmotherm binary file
path_to_cosmotherm = "/Users/boris/Work/SPECTRAL-DISTORTIONS/cosmotherm.rel_corr"

#the path where you want to save the output from cosmotherm
path_to_ct_spectra_results =  "/Users/boris/Work/SPECTRAL-DISTORTIONS/pi_spec/pispec/ct_spectra"
subprocess.call(['mkdir','-p',path_to_ct_spectra_results])
