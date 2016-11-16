from pylada.vasp.relax import Relax
from copy import deepcopy
from pylada.crystal.structure import Structure
from pylada.crystal.atom import Atom
from pylada.vasp import MassExtract,Extract
from pylada.misc import RelativePath
from pylada.vasp.extract import Extract
from pylada.error import ExternalRunFailed
from pylada.vasp import Vasp
import pylada
import os

class CustomChain(object):
    def __init__(self, workflows, names=None, vaspobj=None):
        '''
        Runs a series of workflows
        Args:
            workflow (list): list of lists.  Each sublist should be a series of functions that take, and statically modify, a pylada.vasp.relax.Relax object.
            names (str list):  Titles of runs.  len(workflow) must equal len(names).  Defaults to integers.
            vaspobj: TODO
        '''
        if not names:
            names = [ str(i) for i in range(len(workflows)) ]
        else:
            assert len(names) == len(workflows)
        self.names=names
        self.vasp=vaspobj
        self.workflows = workflows
        self.total_steps = len(workflows)
        self.current_step = 0
        return

    def Extract(self, jobdir):
        extract = MassExtract(jobdir)
        success={}
        for name in self.names:
            #import pdb; pdb.set_trace()
            success[name]=Extract(jobdir+'/'+name).success
        success=all(success.values())
        extract.success=success
        return extract

    def run_calculation(self, name, workflow, vasp, structure, outdir, kwargs):
        structure_ = structure.copy()
        outdir = os.getcwd() if outdir is None else RelativePath(outdir).path
        for modification in workflow:
            vasp = modification(vasp, structure)    
        ## if this calculation has not been done run it
        params = deepcopy(kwargs)
        fulldir = os.path.join(outdir, name)
        output = vasp(structure_, outdir=fulldir, **params)
        if not output.success:
            raise ExternalRunFailed("VASP calculation did not succeed.")

    # Creating the workflow
    def __call__(self, structure, outdir=None, **kwargs ):

        # make this function stateless.
        structure_ = structure.copy()
        for i in range(len(self.workflows)):
            vasp = Relax(copy=deepcopy(self.vasp))
            name = self.names[i]
            workflow = self.workflows[i]
            self.run_calculation(name, workflow, vasp, structure, outdir, kwargs)
        return self.Extract(fulldir)

def load_default_vasp(vasp,structure=None):
    vasp.has_nlep = False
    vasp.maxiter = 100
    if structure==None:
        vasp.kpoints="\n0\nAuto\n14"
    else:
        kpoints_density = 4000
        vasp.kpoints=pylada.gen_kpts(structure,kpoints_density)
    vasp.first_trial = { "kpoints": "\n0\nAuto\n12", "encut": 400.0 }
    vasp.program = '$VASP_KPTS'   # default vasp program
    vasp.npar       = int(structure.__len__() / 24 / 1.25)
    vasp.prec       = "accurate"
    vasp.ediff      = 1.0e-6        # total, not per atom
    vasp.ediffg     = -0.02
    vasp.addgrid    = True
    vasp.ismear     = 0
    vasp.sigma      = 0.01
    vasp.lmaxmix = 4
    vasp.convergence = 1.0e-5
    vasp.algo = "Normal"
    return(vasp)

def all_output(vasp, structure=None):
    vasp.lvtot=True                       #Write LOCPOT
    vasp.add_keyword('laechg','T')   #Print AECCAR for Bader analysis
    vasp.lorbit=11                       #Print PROCAR and DOSCAR
    return vasp

##############
# ELECTRONIC #
##############

def set_algo_normal(vasp: Vasp, structure=None):
    '''

    :param vasp: Vasp
    :param structure:
    :return:
    '''
    vasp.algo = "Normal"
    return vasp

def set_algo_all(vasp: Vasp, structure=None):
    vasp.algo = "All"
    return vasp

def set_algo_fast(vasp: Vasp, structure=None):
    vasp.algo = "Fast"
    return vasp

def hse06(vasp: Vasp, structure=None):
    vasp.istart = 1
    vasp.icharg = 1
    vasp.nelm = 1000
    vasp.add_keyword('lhfcalc', True)
    vasp.precfock = 'Fast'
    vasp.add_keyword('hfscreen', 0.2)
    vasp.algo = 'All'
    return vasp

def set_dos(vasp: Vasp, structure=None):
    vasp.algo = 'None'
    return vasp

def tetrahedron(vasp: Vasp, structure=None):
    vasp.ismear = -5
    return vasp

##########
# IONIC  #
##########

def set_iopt_7(vasp: Vasp, structure=None):
    vasp.ibrion = 3
    vasp.add_keyword('potim', 0)
    vasp.add_keyword('iopt', 7)
    return vasp

def single_point(vasp: Vasp, structure=None):
    vasp.ibrion = -1
    vasp.nsw = 0
    return Vasp

######################
# CONVERGENCE LEVELS #
######################

def awful_converge(vasp: Vasp, structure=None):
    # Start
    vasp.istart = 0
    vasp.icharg = 2
    # Electronic
    vasp.prec = "Normal"
    vasp.nelm = 60
    vasp.nelmin = 8
    vasp.ediff = 1e-4
    # Ionic
    vasp.edffg = -0.2
    # Output
    vasp.lwave = False
    vasp.lcharg = False
    return vasp

def rough_converge(vasp: Vasp, structure=None):
    # Start
    vasp.istart = 0
    vasp.icharg = 2
    # Electronic
    vasp.prec = "Accurate"
    vasp.nelm = 60
    vasp.nelmin = 8
    vasp.ediff = 1e-4
    # Ionic
    vasp.edffg = -0.05
    # Output
    vasp.lwave = False
    vasp.lcharg = False
    return vasp

def get_eigen(vasp: Vasp, structure=None):
    # Start
    vasp.istart = 0
    vasp.icharg = 2
    # Electronic
    vasp.prec = "Accurate"
    vasp.nelm = 200
    vasp.ediff = 1e-6
    # Ionic
    vasp.ediffg = -0.03
    # Output
    vasp.lwave = True
    vasp.lcharg = True
    return vasp

def full_converge(vasp: Vasp, structure=None):
    # Start
    vasp.istart = 1
    vasp.icharg = 1
    # Electronic
    vasp.prec = "Accurate"
    vasp.nelm = 200
    vasp.ediff = 1e-7
    # Ionic
    vasp.ediffg = -0.02
    # Output
    vasp.lwave = True
    vasp.lcharg = True
    return vasp


###########
# KPOINTS #
###########

def set_gamma(vasp: Vasp, structure=None):
    x=1; y=1 ; z=1
    packing = 'Gamma'
    vasp.kpoints = "Gamma_Mesh\n0\n{}\n{} {} {}".format(packing, x, y, z)
    if packing[0].upper() == 'G':
        vasp.add_keyword('auto_gamma', 'True')
    return vasp

def set_222(vasp: Vasp, structure=None):
    x=2; y=2 ; z=2
    packing = 'Gamma'
    vasp.kpoints = "Gamma_Mesh\n0\n{}\n{} {} {}".format(packing, x, y, z)
    if packing[0].upper() == 'G':
        vasp.add_keyword('auto_gamma', 'True')
    return vasp

def set_221(vasp: Vasp, structure=None):
    x=2; y=2 ; z=1
    packing = 'Gamma'
    vasp.kpoints = "Gamma_Mesh\n0\n{}\n{} {} {}".format(packing, x, y, z)
    if packing[0].upper() == 'G':
        vasp.add_keyword('auto_gamma', 'True')
    return vasp