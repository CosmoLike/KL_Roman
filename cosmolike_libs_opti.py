import sys
import emcee
import ctypes
import os
import numpy as np
import mpi4py

# from mpp_blinding import blind_parameters
# from mpp_blinding import seed as blinding_seed

dirname = os.path.split(__file__)[0]
lib_name = os.path.join(dirname, "like_fourier.so")
lib=ctypes.cdll.LoadLibrary(lib_name)
double = ctypes.c_double

Double10 = double*10

initcosmo=lib.init_cosmo_runmode
initcosmo.argtypes=[ctypes.c_char_p]

initcosmo_DEu95CPL=lib.init_cosmo_runmode_DEu95CPL
initcosmo_DEu95CPL.argtypes=[ctypes.c_char_p]

initcosmo_DEl95CPL=lib.init_cosmo_runmode_DEl95CPL
initcosmo_DEl95CPL.argtypes=[ctypes.c_char_p]

initbins=lib.init_binning_fourier
initbins.argtypes=[ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_int]

initsurvey=lib.init_survey
initsurvey.argtypes=[ctypes.c_char_p]

initgalaxies=lib.init_galaxies
initgalaxies.argtypes=[ctypes.c_char_p,ctypes.c_char_p,ctypes.c_char_p,ctypes.c_char_p,ctypes.c_char_p]

initclusters=lib.init_clusters
initclusters.argtypes=[]

initia=lib.init_IA
initia.argtypes=[ctypes.c_char_p,ctypes.c_char_p]

initpriors=lib.init_priors
initpriors.argtypes=[ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p]

initpriors_KL=lib.init_priors_KL
initpriors_KL.argtypes=[ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p]

initpriors_IA_bary=lib.init_priors_IA_bary
initpriors_IA_bary.argtypes=[ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, 
    ctypes.c_bool, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, 
    ctypes.c_bool, ctypes.c_double, ctypes.c_double, ctypes.c_double]

initpriors_IA=lib.init_priors_IA
initpriors_IA.argtypes=[ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]

initprobes=lib.init_probes
initprobes.argtypes=[ctypes.c_char_p]

initdatainv=lib.init_data_inv
initdatainv.argtypes=[ctypes.c_char_p,ctypes.c_char_p]

initdatainvbary=lib.init_data_inv_bary
initdatainvbary.argtypes=[ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p]

get_N_tomo_shear = lib.get_N_tomo_shear
get_N_tomo_shear.argtypes = []
get_N_tomo_shear.restype = ctypes.c_int

get_N_tomo_clustering = lib.get_N_tomo_clustering
get_N_tomo_clustering.argtypes = []
get_N_tomo_clustering.restype = ctypes.c_int

get_N_ggl = lib.get_N_ggl
get_N_ggl.argtypes = []
get_N_ggl.restype = ctypes.c_int

get_N_ell = lib.get_N_ell
get_N_ell.argtypes = []
get_N_ell.restype = ctypes.c_int


# lib.initialize_all_wrapper.restype = ctypes.c_int
# lib.initialize_all_wrapper.argtypes = [
#     ctypes.c_char_p,       # const char * base_dir, 
#     ctypes.c_bool,         # bool auborg_prior, 
#     ctypes.c_bool,         # bool photo_bao_prior, 
#     ctypes.c_bool,         # bool ia_datavector, 
#     ctypes.c_char_p,       # char * ia_model, // must be either "none", "NLA_HF"
#     ctypes.c_char_p,       # char * ia_luminosity_function // should be either "GAMA", "DEEP2"
#     ctypes.c_bool,         # bool modify_shear_priors,
#     Double10,              # double shear_m_mean[10],
#     Double10,              # double shear_m_var[10],
#     ctypes.c_bool,         # bool modify_photoz_priors,
#     Double10,              # double photoz_source_bias_mean[10],
#     double,                # double photoz_source_sigma_mean,
#     Double10,              # double photoz_source_bias_var[10],
#     double,                # double photoz_source_sigma_var,
#     ctypes.c_char_p,       # source_filename
#     ctypes.c_char_p,       # lens_filename
# ]


# initialize_all_wrapper=lib.initialize_all_wrapper

class IterableStruct(ctypes.Structure):
    def names(self):
        out = []
        for name, obj, length in self.iter_parameters():
            if length==0:
                out.append(name)
            else:
                for i in xrange(length):
                    out.append(name + "_" + str(i))
        return out


    def iter_parameters(self):
        for name,ptype in self._fields_:
            obj = getattr(self, name)
            if hasattr(ptype, "_length_"):
                yield name, obj, ptype._length_
            else:
                yield name, obj, 0

    def iter_parameters_filter(self, used):
        for (name, obj, length) in self.iter_parameters():
            if name in used:
                yield name, obj, 0


    def convert_to_vector(self):
        p = []
        for name, obj, length in self.iter_parameters():
            if length==0:
                p.append(obj)
            else:
                for i in xrange(length):
                    p.append(obj[i])
        return p

    def convert_to_vector_filter(self, used):
        p = []
        for name, obj, length in self.iter_parameters():
            if length==0:
                if name in used:
                    p.append(obj)
            else:
                for i in xrange(length):
                    if name+'_'+str(i) in used:
                        p.append(obj[i])
        return p



    def read_from_cosmosis(self, block):
        for name,ptype in self._fields_:
            obj = getattr(self, name)
            if hasattr(ptype, "_length_"):
                for i in xrange(ptype._length_):
                    obj[i] = block[self.section_name, name+"_"+str(i)]
            else:
                setattr(self, name, block[self.section_name, name])



    def print_struct(self):
        for name,ptype in self._fields_:
            obj = getattr(self, name)
            if hasattr(ptype, "_length_"):
                for i in xrange(ptype._length_):
                    print "%s[%d] = %f" % (name, i, obj[i])
            else:
                print "%s = %f" % (name, obj)


    def number_of_doubles(self):
        n=0
        for name, ptype in self._fields_:
            if hasattr(ptype, "_length_"):
                n += ptype._length_
            else:
                n += 1
        return n

    def set_from_vector(self, p):
        i=0
        j=0
        while i<len(p):
            name,ptype = self._fields_[j]
            j+=1
            if ptype == double:
                setattr(self, name, p[i])
                i+=1
            else:
                x = getattr(self, name)
                assert x._type_==double
                for k in xrange(x._length_):
                    x[k] = p[i]
                    i+=1


class InputCosmologyParams(IterableStruct):
    section_name = "cosmological_parameters"
    _fields_ = [
        ("omega_m", double),
        ("sigma_8", double),
        ("n_s", double),
        ("w0", double),
        ("wa", double),
        ("omega_b", double),
        ("h0", double),
        ("MGSigma", double),
        ("MGmu", double),
    ]

    @classmethod
    def fiducial(cls):
        c = cls()
        c.omega_m = 0.3156
        c.sigma_8 = 0.831
        c.n_s = 0.9645
        c.w0 = -1.0
        c.wa = 0.0
        c.omega_b = 0.0491685
        c.h0 = 0.6727
        c.MGSigma = 0.0
        c.MGmu = 0.0
        return c
    
    @classmethod
    def DEu95CPL(cls):
        c = cls()
        c.omega_m = 0.3156
        c.sigma_8 = 0.831
        c.n_s = 0.9645
        c.w0 = -1.249
        c.wa = 0.59
        c.omega_b = 0.0491685
        c.h0 = 0.6727
        c.MGSigma = 0.0
        c.MGmu = 0.0
        return c

    @classmethod
    def DEl95CPL(cls):
        c = cls()
        c.omega_m = 0.3156
        c.sigma_8 = 0.831
        c.n_s = 0.9645
        c.w0 = -0.289
        c.wa = -2.21
        c.omega_b = 0.0491685
        c.h0 = 0.6727
        c.MGSigma = 0.0
        c.MGmu = 0.0
        return c

    # @classmethod
    # def fiducial_sigma(cls):
    #     c = cls()
    #     c.omega_m = 0.01
    #     c.sigma_8 = 0.01
    #     c.n_s = 0.01
    #     c.w0 = .02
    #     c.wa = 0.02
    #     c.omega_b = 0.001
    #     c.h0 = 0.01
    #     c.MGSigma = 0.1
    #     c.MGmu = 0.1     
    #     return c

    @classmethod
    def fiducial_sigma(cls):
        c = cls()
        c.omega_m = 0.1
        c.sigma_8 = 0.1
        c.n_s = 0.05
        c.w0 = 0.2
        c.wa = 0.2
        c.omega_b = 0.01
        c.h0 = 0.1
        c.MGSigma = 0.1
        c.MGmu = 0.2     
        return c

class InputCosmologyLowzParams(IterableStruct):
    section_name = "cosmological_parameters_lowz"
    _fields_ = [
        ("sigma_8_lowz", double),
        ("z_low", double),
    ]

    @classmethod
    def fiducial(cls):
        c = cls()
        c.sigma_8_lowz = 0.831
        c.z_low = 1.0
        return c

    @classmethod
    def fiducial_sigma(cls):
        c = cls()
        c.sigma_8_lowz = 0.0001
        c.z_low = 0.02
        return c

class InputNuisanceParams(IterableStruct):
    section_name = "nuisance_parameters"
    _fields_ = [
        ("bias", double*10),
        ("source_z_bias", double*10),
        ("source_z_s", double),
        ("lens_z_bias", double*10),
        ("lens_z_s", double),
        ("shear_m", double*10),
        ("A_ia", double),
        ("beta_ia", double),
        ("eta_ia", double),
        ("eta_ia_highz", double),
        ("lf", double*6),
        ("m_lambda", double*6),
        ("cluster_c", double*4),
        ("bary", double*3)
    ]

    @classmethod
    def fiducial(cls):
        c = cls()
        c.bias[:] = [1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2]
        c.source_z_bias[:] = np.repeat(0.0, 10)
        c.source_z_s = 0.05
        c.lens_z_bias[:] = np.repeat(0.0, 10)
        c.lens_z_s = 0.03
        c.shear_m[:] = np.repeat(0.0, 10)
        c.A_ia = 5.92# 5.95 from Eifler et. al. 2020
        c.beta_ia = 1.1# 1.1
        c.eta_ia = -0.47# 0.49
        c.eta_ia_highz = 0.0# 0.0
        c.lf[:] = np.repeat(0.0, 6)
        c.m_lambda[:] = [3.207, 0.993, 0.0, 0.456, 0.0, 0.0]
        c.cluster_c[:] = [0., 0., 0., 0.]
        c.bary[:] = [0., 0., 0.]
        return c

    @classmethod
    def fiducial_KL(cls):
        c = cls()
        c.bias[:] = [1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75]
        c.source_z_bias[:] = np.repeat(0.0, 10)
        c.source_z_s = 0.0005
        c.lens_z_bias[:] = np.repeat(0.0, 10)
        c.lens_z_s = 0.0005
        c.shear_m[:] = np.repeat(0.0, 10)
        c.A_ia = 5.92
        c.beta_ia = 1.1
        c.eta_ia = -0.47
        c.eta_ia_highz = 0.0
        c.lf[:] = np.repeat(0.0, 6)
        c.m_lambda[:] = [3.207, 0.993, 0.0, 0.456, 0.0, 0.0]
        c.cluster_c[:] = [0., 0., 0., 0.]
        c.bary[:] = [0., 0., 0.]
        return c

    
    @classmethod
    def fiducial_sigma(cls):
        c = cls()
        c.bias[:] = np.repeat(0.15, 10)
        c.source_z_bias[:] = np.repeat(0.005, 10)
        c.source_z_s = 0.002
        c.lens_z_bias[:] = np.repeat(0.005, 10)
        c.lens_z_s = 0.002
        c.shear_m[:] = np.repeat(0.005, 10)
        c.A_ia = 0.05
        c.beta_ia = 0.01
        c.eta_ia = 0.01
        c.eta_ia_highz = 0.01
        c.lf[:] = np.repeat(0.005, 6)
        c.m_lambda[:] = [0.045, 0.045, 0.3, 0.045, 0.03, 0.1]
        c.cluster_c[:] = [1., 1., 1., 1.]
        c.bary[:] = [3., 1., .15]
        return c        

    @classmethod
    def fiducial_sigma_KL(cls):
        c = cls()
        c.bias[:] = np.repeat(0.15, 10)
        c.source_z_bias[:] = np.repeat(0.00025, 10)
        c.source_z_s = 0.0001
        c.lens_z_bias[:] = np.repeat(0.00025, 10)
        c.lens_z_s = 0.0001
        c.shear_m[:] = np.repeat(0.005, 10)
        c.A_ia = 0.05
        c.beta_ia = 0.01
        c.eta_ia = 0.01
        c.eta_ia_highz = 0.01
        c.lf[:] = np.repeat(0.005, 6)
        c.m_lambda[:] = [0.045, 0.045, 0.3, 0.045, 0.03, 0.1]
        c.cluster_c[:] = [1., 1., 1., 1.]
        c.bary[:] = [3., 1., .15]
        return c        

##########

class InputNuisanceParamsGRS(IterableStruct):
    section_name = "nuisance_parameters_GRS"
    _fields_ = [
        ("grsbias", double*7),
        ("grssigmap", double*7),
        ("grssigmaz", double),
        ("grspshot", double),
        ("grskstar", double),
    ]

    @classmethod
    def fiducial(cls):
        c = cls()
        c.grsbias[:] = [1.538026692020565,1.862707210288686,2.213131761595241,2.617023657038295,
                        2.975011712138650,3.376705680190931,3.725882076395691]
        c.grssigmap[:] = np.repeat(290.,7)
        c.grssigmaz = 0.001
        c.grspshot = 0.0
        c.grskstar = 0.24
        return c

    @classmethod
    def fiducial_sigma(cls):
        c = cls()
        c.grsbias[:] = np.repeat(0.15, 7)
        c.grssigmap[:] = np.repeat(20.0, 7)
        c.grssigmaz = 0.0002 # fid is 0.001 and can't be neg
        c.grspshot = 0.001 #fid is zero
        c.grskstar = 0.05 # fid is 0.24
        return c

class LikelihoodFunctionWrapper(object):
    def __init__(self, varied_parameters, KL=False):
        self.varied_parameters = varied_parameters
        if KL:
            self.KL_flag = True
        else:
            self.KL_flag = False


    def fill_varied(self, icp, inp, inpgrs, icplow, x):
        assert len(x) == len(self.varied_parameters), "Wrong number of parameters"
        i = 0
        #print "filling %d params"%(len(x))
        # THE SEQUENCE IS IMPORTANT! MATCH IT WITH YOUR STARING POINT SEQUENCE AND VARIAED
        # PARAMETERS!
        for s in [icp, icplow, inp, inpgrs]:
            for name, obj, length in s.iter_parameters():
                if length==0:
                    if name in self.varied_parameters:
                        setattr(s, name, x[i])
                        #print "Found %s! setting %s to %e [%d/%d]"%(name, name, x[i], i+1, len(x))
                        i+=1
                    #else:
                        #print "%s not found, next one..."%name
                else:
                    for j in xrange(length):
                        name_i = name + "_" + str(j)
                        if name_i in self.varied_parameters:
                            #print "Found %s! setting %s to %e [%d/%d]"%(name_i, name_i, x[i], i+1, len(x))
                            obj[j] = x[i]
                            i+=1
                        #else:
                            #print "%s not found, next one..."%name_i

    def __call__(self, x):
        icp = InputCosmologyParams.fiducial()
        if self.KL_flag:
            inp = InputNuisanceParams.fiducial_KL()
        else:
            inp = InputNuisanceParams.fiducial()
        inpgrs = InputNuisanceParamsGRS.fiducial()
        icplow = InputCosmologyLowzParams.fiducial()
        # fill the param-struct with input values and default values
        self.fill_varied(icp, inp, inpgrs, icplow, x)
        #icp.print_struct()
        #inp.print_struct()
        #inpgrs.print_struct()
        #icplow.print_struct()
        #print
        like = lib.log_like_wrapper(icp, inp, inpgrs, icplow)
        #print "like before" , like
        if like < -1.0e+14:
            return -np.inf
        return like


lib.log_like_wrapper.argtypes = [InputCosmologyParams, InputNuisanceParams, InputNuisanceParamsGRS, InputCosmologyLowzParams]
lib.log_like_wrapper.restype = double
log_like_wrapper = lib.log_like_wrapper


def sample_cosmology_only_SN_SL(MG = False):
    if MG:
        varied_parameters = InputCosmologyParams().names()
    else:
        varied_parameters = ['omega_m']
        varied_parameters.append('w0')
        varied_parameters.append('wa')
        varied_parameters.append('h0')

    return varied_parameters

def sample_LCDM_only(MG = False):
    if MG:
        varied_parameters = InputCosmologyParams().names()
    else:
        varied_parameters = ['omega_m']
        varied_parameters.append('sigma_8')
        varied_parameters.append('n_s')
        varied_parameters.append('omega_b')
        varied_parameters.append('h0')

    return varied_parameters

def sample_LCDM_sigma8_split(MG = False):
    if MG:
        varied_parameters = InputCosmologyParams().names()
    else:
        varied_parameters = ['omega_m']
        varied_parameters.append('sigma_8')
        varied_parameters.append('n_s')
        varied_parameters.append('omega_b')
        varied_parameters.append('h0')
        varied_parameters.append('sigma_8_lowz')
    return varied_parameters

def sample_sigma8_split():
    varied_parameters = ['sigma_8', 'sigma_8_lowz']
    return varied_parameters

def sample_cosmology_only(MG = False):
    if MG:
        varied_parameters = InputCosmologyParams().names()
    else:
        varied_parameters = ['omega_m']
        varied_parameters.append('sigma_8')
        varied_parameters.append('n_s')
        varied_parameters.append('w0')
        varied_parameters.append('wa')
        varied_parameters.append('omega_b')
        varied_parameters.append('h0')
    return varied_parameters

def sample_cosmology_shear_nuisance_IA(tomo_N_shear,MG = False, DE=True):
    if DE:
        varied_parameters = sample_cosmology_only(MG)
    else:
        varied_parameters = sample_LCDM_only(MG)
    varied_parameters += ['source_z_bias_%d'%i for i in xrange(tomo_N_shear)]
    varied_parameters.append('source_z_s')
    varied_parameters += ['shear_m_%d'%i for i in xrange(tomo_N_shear)]
    varied_parameters.append('A_ia')
    varied_parameters.append('beta_ia')
    varied_parameters.append('eta_ia')
    varied_parameters.append('eta_ia_highz')
    return varied_parameters

def sample_cosmology_shear_nuisance_IA_bary(tomo_N_shear,MG = False, NPCs=3, DE=True):
    if DE:
        varied_parameters = sample_cosmology_only(MG)
    else:
        varied_parameters = sample_LCDM_only(MG)
    varied_parameters += ['source_z_bias_%d'%i for i in xrange(tomo_N_shear)]
    varied_parameters.append('source_z_s')
    varied_parameters += ['shear_m_%d'%i for i in xrange(tomo_N_shear)]
    varied_parameters.append('A_ia')
    varied_parameters.append('beta_ia')
    varied_parameters.append('eta_ia')
    varied_parameters.append('eta_ia_highz')
    varied_parameters += ['bary_%d'%i for i in xrange(NPCs)]
    return varied_parameters

def sample_cosmology_shear_nuisance(tomo_N_shear, MG=False, NPCs=0, cosmology='LCDM', source_photo_z=True, shear_calibration=True, IA=True):
    if(cosmology=='LCDM'):
        varied_parameters = sample_LCDM_only(MG)
    elif(cosmology=='wCDM'):
        varied_parameters = sample_LCDM_only(MG)
        varied_parameters.append("w0")
    elif(cosmology=='w0waCDM'):
        varied_parameters = sample_cosmology_only(MG)
    elif(cosmology=='LCDM_split'):
        varied_parameters = sample_LCDM_sigma8_split(MG)
    elif(cosmology=='wCDM_split'):
        varied_parameters = sample_LCDM_sigma8_split(MG)
        varied_parameters.append("w0")
    elif(cosmology=='w0waCDM_split'):
        varied_parameters = sample_LCDM_sigma8_split(MG)
        varied_parameters.append('w0')
        varied_parameters.append('wa')
    elif(cosmology=='s8split_only'):
        varied_parameters = sample_sigma8_split()
    else:
        print "cosmology = %s is not supported!"%cosmology
        exit(-1)
    if source_photo_z:
        varied_parameters+=['source_z_bias_%d'%i for i in xrange(tomo_N_shear)]
        varied_parameters.append('source_z_s')
    if shear_calibration:
        varied_parameters+=['shear_m_%d'%i for i in xrange(tomo_N_shear)]
    if IA:
        varied_parameters.append('A_ia')
        varied_parameters.append('beta_ia')
        varied_parameters.append('eta_ia')
        varied_parameters.append('eta_ia_highz')
    if NPCs>0:
        varied_parameters += ['bary_%d'%i for i in xrange(NPCs)]
    return varied_parameters

def sample_cosmology_clustering_nuisance(tomo_N_lens,MG = False):
    varied_parameters = sample_cosmology_only(MG)
    varied_parameters += ['bias_%d'%i for i in xrange(tomo_N_lens)]
    varied_parameters += ['lens_z_bias_%d'%i for i in xrange(tomo_N_lens)]
    varied_parameters.append('lens_z_s')
    return varied_parameters

def sample_cosmology_2pt_nuisance(tomo_N_shear,tomo_N_lens,MG = False):
    varied_parameters = sample_cosmology_only(MG)
    varied_parameters += ['bias_%d'%i for i in xrange(tomo_N_lens)]
    varied_parameters += ['source_z_bias_%d'%i for i in xrange(tomo_N_shear)]
    varied_parameters.append('source_z_s')
    varied_parameters += ['lens_z_bias_%d'%i for i in xrange(tomo_N_lens)]
    varied_parameters.append('lens_z_s')
    varied_parameters += ['shear_m_%d'%i for i in xrange(tomo_N_shear)]
    return varied_parameters

def sample_cosmology_2pt_shear_nuisance(tomo_N_shear,tomo_N_lens,MG = False):
    varied_parameters = sample_cosmology_only(MG)
    varied_parameters += ['bias_%d'%i for i in xrange(tomo_N_lens)]
    varied_parameters += ['shear_m_%d'%i for i in xrange(tomo_N_shear)]
    return varied_parameters

def sample_cosmology_2pt_photo_nuisance(tomo_N_shear,tomo_N_lens,MG = False):
    varied_parameters = sample_cosmology_only(MG)
    varied_parameters += ['bias_%d'%i for i in xrange(tomo_N_lens)]
    varied_parameters += ['source_z_bias_%d'%i for i in xrange(tomo_N_shear)]
    varied_parameters.append('source_z_s')
    varied_parameters += ['lens_z_bias_%d'%i for i in xrange(tomo_N_lens)]
    varied_parameters.append('lens_z_s')
    return varied_parameters


def sample_cosmology_2pt_cluster_nuisance(tomo_N_shear,tomo_N_lens,MG = False):
    if MG:
        print "sample_cosmology_2pt_cluster_nuisance: MG = True not yet supported for clusters"
        os.exit()
    varied_parameters = sample_cosmology_only(MG)
    varied_parameters += ['bias_%d'%i for i in xrange(tomo_N_lens)]
    varied_parameters += ['source_z_bias_%d'%i for i in xrange(tomo_N_shear)]
    varied_parameters.append('source_z_s')
    varied_parameters += ['lens_z_bias_%d'%i for i in xrange(tomo_N_lens)]
    varied_parameters.append('lens_z_s')
    varied_parameters += ['shear_m_%d'%i for i in xrange(tomo_N_shear)]
    varied_parameters += ['m_lambda_%d'%i for i in xrange(6)]
    return varied_parameters

def sample_cosmology_clusterN_clusterWL_nuisance(tomo_N_shear,MG = False):
    if MG:
        print "sample_cosmology_clusterN_clusterWL_nuisance: MG = True not yet supported for clusters"
        os.exit()
    varied_parameters = sample_cosmology_only(MG)
    varied_parameters += ['source_z_bias_%d'%i for i in xrange(tomo_N_shear)]
    varied_parameters.append('source_z_s')
    varied_parameters += ['shear_m_%d'%i for i in xrange(tomo_N_shear)]
    varied_parameters += ['m_lambda_%d'%i for i in xrange(6)]
    return varied_parameters

def sample_cosmology_clusterN_clusterWL_nuisance_fixMOR(tomo_N_shear,MG = False):
    if MG:
        print "sample_cosmology_clusterN_clusterWL_nuisance: MG = True not yet supported for clusters"
        os.exit()
    varied_parameters = sample_cosmology_only(MG)
    varied_parameters += ['source_z_bias_%d'%i for i in xrange(tomo_N_shear)]
    varied_parameters.append('source_z_s')
    varied_parameters += ['shear_m_%d'%i for i in xrange(tomo_N_shear)]
    return varied_parameters


def sample_cosmology_2pt_cluster_nuisance_GRS(tomo_N_shear,tomo_N_lens,MG = False):
    if MG:
        print "sample_cosmology_2pt_cluster_nuisance: MG = True not yet supported for clusters"
        os.exit()
    varied_parameters = sample_cosmology_only(MG)
    varied_parameters += ['bias_%d'%i for i in xrange(tomo_N_lens)]
    varied_parameters += ['source_z_bias_%d'%i for i in xrange(tomo_N_shear)]
    varied_parameters.append('source_z_s')
    varied_parameters += ['lens_z_bias_%d'%i for i in xrange(tomo_N_lens)]
    varied_parameters.append('lens_z_s')
    varied_parameters += ['shear_m_%d'%i for i in xrange(tomo_N_shear)]
    varied_parameters += ['m_lambda_%d'%i for i in xrange(6)]
    varied_parameters += ['grsbias_%d'%i for i in xrange(7)]
    varied_parameters += ['grssigmap_%d'%i for i in xrange(7)]
    varied_parameters.append("grssigmaz")
    varied_parameters.append("grskstar")
    return varied_parameters

def sample_cosmology_grs_nuisance(MG = False):
    varied_parameters = sample_cosmology(MG)
    varied_parameters += ['grsbias_%d'%i for i in xrange(7)]
    varied_parameters += ['grssigmap_%d'%i for i in xrange(7)]
    varied_parameters.append("grssigmaz")
    varied_parameters.append("grskstar")

    return varied_parameters

def sample_cosmology_SN_WFIRST():
    varied_parameters = ['omega_m']
    varied_parameters.append('w0')
    varied_parameters.append('wa')
    varied_parameters.append('h0')

    return varied_parameters

def sample_main(varied_parameters, iterations, nwalker, nthreads, filename, blind=False, pool=None, KL=False):
    print varied_parameters

    ### Choose your cosmology: Fiducial or DEu/l95CPL?
    #likelihood = LikelihoodFunctionWrapper(varied_parameters)
    starting_point = InputCosmologyParams.fiducial().convert_to_vector_filter(varied_parameters)
    #starting_point += InputNuisanceParams().fiducial().convert_to_vector_filter(varied_parameters)

    std = InputCosmologyParams.fiducial_sigma().convert_to_vector_filter(varied_parameters)
    #std += InputNuisanceParams().fiducial_sigma().convert_to_vector_filter(varied_parameters)

    ### Choose your nuisance: weak lensing or kinematic lensing?
    if KL:
        starting_point += InputCosmologyLowzParams().fiducial().convert_to_vector_filter(varied_parameters)
        starting_point += InputNuisanceParams().fiducial_KL().convert_to_vector_filter(varied_parameters)
        starting_point += InputNuisanceParamsGRS().fiducial().convert_to_vector_filter(varied_parameters)
        std += InputCosmologyLowzParams().fiducial_sigma().convert_to_vector_filter(varied_parameters)
        std += InputNuisanceParams().fiducial_sigma_KL().convert_to_vector_filter(varied_parameters)
        std += InputNuisanceParamsGRS().fiducial_sigma().convert_to_vector_filter(varied_parameters)
        likelihood = LikelihoodFunctionWrapper(varied_parameters,KL=True)
    else:
        starting_point += InputCosmologyLowzParams().fiducial().convert_to_vector_filter(varied_parameters)
        starting_point += InputNuisanceParams().fiducial().convert_to_vector_filter(varied_parameters)
        starting_point += InputNuisanceParamsGRS().fiducial().convert_to_vector_filter(varied_parameters)
        std += InputCosmologyLowzParams().fiducial_sigma().convert_to_vector_filter(varied_parameters)
        std += InputNuisanceParams().fiducial_sigma().convert_to_vector_filter(varied_parameters)
        std += InputNuisanceParamsGRS().fiducial_sigma().convert_to_vector_filter(varied_parameters)
        likelihood = LikelihoodFunctionWrapper(varied_parameters)

    p0 = emcee.utils.sample_ball(starting_point, std/10000, size=nwalker)

    ndim = len(starting_point)
    print "ndim = ", ndim
    print "start = ", starting_point
    print "std = ", std
    print "p0 = ", p0
    #exit(0)
    #print likelihood(starting_point)
    #exit(0)
    if pool is not None:
        if not pool.is_master():
            pool.wait()
            sys.exit(0)


    sampler = emcee.EnsembleSampler(nwalker, ndim, likelihood,threads=nthreads,pool=pool)

#    sampler = emcee.EnsembleSampler(nwalker, ndim, likelihood, pool=pool)

    f = open(filename, 'w')
    if f:
        print "Writing output file..."
    else:
        print "Fail to open %s"%filename
    #write header here
    f.write('# ' + '    '.join(varied_parameters)+" log_like\n")
    f.write('#blind=%s\n'%blind)
    if blind:
        f.write('#blinding_seed=%d\n'%blinding_seed)

    for (p, loglike, state) in sampler.sample(p0,iterations=iterations):
        for row,logl in zip(p,loglike):
            if blind:
                row = blind_parameters(varied_parameters, row)
            p_text = '  '.join(str(r) for r in row)
            #print p_text, row
            f.write('%s %e\n' % (p_text,logl))
        f.flush()
    f.close()
    
    pool.close()

    # for (p, loglike, state) in sampler.sample(p0,iterations=iterations):
    #     for row in p:
    #         if blind:
    #             row = blind_parameters(varied_parameters, row)
    #         p_text = '  '.join(str(r) for r in row)
    #         print ('%s %e\n' % (p_text,loglike))
    #         f.write('%s %e\n' % (p_text,loglike))
    #     f.flush()
    # f.close()

