from .utils import *
from .standard_mu_and_y_distortions import *
from .Recfast_wrapper import *

class Xe_PCA_EigenModes:
    def __init__(self):
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

        E1 = np.loadtxt(path_to_ct_database+'../PCA_modes/Modes/mode_N121_so_planck_1.dat')
        E2 = np.loadtxt(path_to_ct_database+'../PCA_modes/Modes/mode_N121_so_planck_2.dat')
        E3 = np.loadtxt(path_to_ct_database+'../PCA_modes/Modes/mode_N121_so_planck_3.dat')
        self.Xe_PCA_EigenModes['E1']['z'] = E1[:,0]
        self.Xe_PCA_EigenModes['E2']['z'] = E2[:,0]
        self.Xe_PCA_EigenModes['E3']['z'] = E3[:,0]
        self.Xe_PCA_EigenModes['E1']['values'] = E1[:,1]
        self.Xe_PCA_EigenModes['E2']['values'] = E2[:,1]
        self.Xe_PCA_EigenModes['E3']['values'] = E3[:,1]
        self.Xe_PCA_EigenModes['E1']['sigma1'] = 0.12
        self.Xe_PCA_EigenModes['E2']['sigma2'] = 0.19
        self.Xe_PCA_EigenModes['E3']['sigma3'] = 0.35
