Python package to study spectral distortions of the cosmic microwave background radiation.

Mainly based on codes developed by Jens Chluba and the Manchester spectral distortion group (Sandeep Acharya, Boris Bolliet, Luke Hart, Charis Kaur, Tom Kite, Elizabeth Lee, Mathieu Remazeilles, Francesco Pace, Andrea Ravenni, Aditya Rotti, Abir Sarkar).

User dependent paths have to be set in specdist/config.py.

To learn about CMB spectral distortions:

https://physique.cuso.ch/fileadmin/physique/document/2014_Chluba_notes.pdf


https://arxiv.org/pdf/1909.01593.pdf


Cosmotherm comments:

* redshift dependent functions requiring PDE solver in output_fz_functions(), inside Solve_PDE.cpp.

* photon injection source and heating term in cosmotherm.rel_corr/Thermalization_Module/define_PDE_Plugins/define_PDE_photon_injection_decay.cpp

* temperature equation (rho_e) defined in def_drho_dz in define_PDE.cpp
