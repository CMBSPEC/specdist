from .utils import *
from .standard_mu_and_y_distortions import *

class edges:
    edges_nu_in_MHz = 78.
    edges_nu = edges_nu_in_MHz*1e6
    edges_x = hplanck*edges_nu/kb/firas_T0_bf
    edges_I = B_nu_of_T(edges_nu,firas_T0)
