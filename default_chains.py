from pylada.vasp.relax import Relax
from pylada.vasp.functional import Vasp
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
import math


class CustomFunctional(object):
    def __init__(self, base: Vasp, modifications: list):
        self.base = base
        self.modifications = modifications
        return

class CustomChain(object):
    def __init__(self, functionals: list, names=None, vaspobj:Vasp=None, basename=''):
        '''
        Runs a series of workflows
        Args:
            workflow (list): list of lists.  Each sublist should be a series of functions that take, and statically modify, a pylada.vasp.relax.Relax object.
            names (str list):  Titles of runs.  len(workflow) must equal len(names).  Defaults to integers.
            vaspobj: TODO
        '''
        if not names:
            names = [str(i) for i in range(len(functionals))]
        else:
            assert len(names) == len(functionals)
        self.names=names
        self.vasp = vaspobj
        self.functionals = functionals
        self.total_steps = len(functionals)
        self.current_step = 0
        self.basename = basename
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

    def run_calculation(self, name, workflow: CustomFunctional, structure, outdir, previous, **kwargs):
        vasp = workflow.base(copy=deepcopy(self.vasp))
        structure_ = structure.copy()
        outdir = os.getcwd() if outdir is None else RelativePath(outdir).path
        for modification in workflow.modifications:
            vasp = modification(vasp, structure_)
        ## if this calculation has not been done run it
        params = deepcopy(kwargs)
        fulldir = os.path.join(outdir, name)
        output = vasp(structure_, outdir=fulldir, restart=previous ,**params)
        if not output.success:
            raise ExternalRunFailed("VASP calculation did not succeed.")
        return output

    # Creating the workflow
    def __call__(self, structure, outdir=None, names=None, functionals=None, **kwargs ):
        if not names:
            names = self.names
        if not functionals:
            functionals = self.functionals
        # make this function stateless.
        previous =  None
        for i in range(len(functionals)):
            name = names[i]
            workflow = functionals[i]
            print(previous)
            previous = self.run_calculation(name, workflow, structure, outdir, previous, **kwargs)
        fulldir = os.path.join(outdir, name)
        return self.Extract(fulldir)

class SpinCustomChain(CustomChain):
    def __init__(self, functionals: list, nupdown_functionals : list, nupdowns, names=None, vaspobj:Vasp=None, basename=''):
        '''
        Runs a series of workflows
        Args:
            functionals (list): list of lists.  Each sublist should be a series of functions that take, and statically modify, a pylada.vasp.relax.Relax object.
            nupdown_functionals (list) : list of lists.  Each sublist should be a series of functions that take, and statically modify, a pylada.vasp.relax.Relax object.
            nupdowns (int list): nupdowns that will be considered
            vaspobj: TODO
        '''
        self.nupdowns =  nupdowns
        self.nupdown_functionals = nupdown_functionals
        return  super().__init__(functionals, names=names, vaspobj=vaspobj, basename=basename)

    def __call__(self, structure, outdir=None, **kwargs):
        energies = {}
        for nup in self.nupdowns:  # Check energies of various spin configurations
            nupdown_outdir = os.path.join(outdir, str(nup))
            names = [ str(x) for x in range(len(self.nupdown_functionals)) ]
            try: # Load answer from directory if it is present
                energies[nup] = float(self.Extract(os.path.join(nupdown_outdir, names[-1])).energy)
                break
            except:  # if the directory does not have the information, run vasp
                def set_nupdown(vasp: Vasp, structure=None):
                    vasp.nupdown = nup
                    return vasp
                for x in self.nupdown_functionals: # Set nupdown
                    x.modifications.append(set_nupdown)
                super().__call__(structure, outdir=nupdown_outdir, functionals=self.nupdown_functionals, names=names)
                energies[nup] = float(self.Extract(os.path.join(nupdown_outdir, names[-1])).energy)

        nup = min(energies.iterkeys(), key=(lambda key: energies[key]))
        def set_nupdown(vasp: Vasp, structure=None):
            vasp.nupdown = nup
            return vasp
        for x in self.functionals: # Set nupdown
            x.modifications.append(set_nupdown)
        return super().__call__(structure, outdir=outdir, **kwargs)

def load_default_vasp(vasp: Vasp,structure=None):
    vasp.has_nlep = False
    vasp.maxiter = 100
    if structure==None:
        vasp.kpoints="\n0\nAuto\n14"
    else:
        kpoints_density = 4000
        vasp.kpoints=pylada.gen_kpts(structure,kpoints_density)
    vasp.first_trial = { "kpoints": "\n0\nAuto\n12", "encut": 400.0 }
    vasp.program = '$VASP_PYLADA'   # default vasp program
    vasp.npar       = int(os.environ['PBS_NUM_NODES'])
    vasp.prec       = "accurate"
    vasp.ediff      = 1.0e-6        # total, not per atom
    vasp.ediffg     = -0.02
    vasp.addgrid    = True
    vasp.ismear     = 0
    vasp.sigma      = 0.01
    vasp.lmaxmix = 4
    vasp.convergence = 1.0e-5
    vasp.algo = "Normal"
    vasp.lorbit=11
    return(vasp)

def all_output(vasp, structure=None):
    vasp.lvtot=True                       #Write LOCPOT
    vasp.add_keyword('lvhar', True)
    vasp.add_keyword('laechg',True)   #Print AECCAR for Bader analysis
    vasp.lorbit=11                       #Print PROCAR and DOSCAR
    return vasp

def set_kpar_auto(vasp: Vasp, structure=None):
    nodes = int(os.environ['PBS_NUM_NODES'])
    np = int(os.environ['PBS_NP'])
    atoms = len(structure)
    if np / atoms > 1:
        kpar = math.ceil(np/atoms)
        kpoint_str = vasp.kpoints.split('\n')[3]
        kpoints = [ int(x) for x in kpoint_str.split() ]
        num_kpoints = kpoints[0] * kpoints[1] * kpoints[2]
        while kpar > 1:
            if num_kpoints % kpar == 0 and nodes % kpar ==0:
                vasp.add_keyword('kpar', kpar)
                vasp.npar = int(nodes / kpar)
                return vasp
            else:
                kpar = kpar -1

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

def set_algo_conj(vasp: Vasp, structure=None):
    vasp.algo = "Conjugate"
    return vasp

def hse06(vasp: Vasp, structure=None):
    vasp.istart = 1
    vasp.icharg = 1
    vasp.nelmdl = 0
    vasp.nelm = 1000
    vasp.add_keyword('lhfcalc', True)
    vasp.precfock = 'Fast'
    vasp.add_keyword('hfscreen', 0.2)
    vasp.algo = 'All'
    vasp.ismear = 0
    vasp.npar = None
    vasp.ldau = False
    vasp.lwave = True
    vasp.lcharg = True
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
    vasp.add_keyword('iopt', None)
    vasp.ediff = 1e-5
    vasp.add_keyword('lmaxmix', None)
    return vasp

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
    vasp.ediff = 1e-3
    vasp.nelmdl = -15
    # Ionic
    vasp.nsw = 5000
    vasp.ediffg = -0.5
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
    vasp.ediff = 5e-4
    vasp.nelmdl = -15
    # Ionic
    vasp.nsw = 5000
    vasp.ediffg = -0.1
    # Output
    vasp.lwave = False
    vasp.lcharg = False
    return vasp

def get_eigen_nospin(vasp: Vasp, structure=None):
    # Start
    vasp.istart = 0
    vasp.icharg = 2
    vasp.ispin = 1
    # Electronic
    vasp.prec = "Accurate"
    vasp.nelm = 200
    vasp.ediff = 5e-4
    vasp.nelmdl = -15
    # Ionic
    vasp.nsw = 5
    vasp.ediffg = -1
    # Output
    vasp.lwave = True
    vasp.lcharg = True
    return vasp

def get_eigen(vasp: Vasp, structure=None):
    # Start
    vasp.istart = 1
    vasp.icharg = 1
    # Electronic
    vasp.prec = "Accurate"
    vasp.nelm = 200
    vasp.ediff = 1e-5
    vasp.nelmdl = 0
    # Ionic
    vasp.nsw = 5000
    vasp.ediffg = -0.05
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
    vasp.ediff = 1e-5
    vasp.nelmdl = 0
    # Ionic
    vasp.nsw = 5000
    vasp.ediffg = -0.03
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
    vasp.ismear = 0
    vasp.add_keyword('kpar', 1)
    return vasp

def set_222(vasp: Vasp, structure=None):
    x=2; y=2 ; z=2
    packing = 'Gamma'
    vasp.kpoints = "Gamma_Mesh\n0\n{}\n{} {} {}".format(packing, x, y, z)
    return vasp

def set_221(vasp: Vasp, structure=None):
    x=2; y=2 ; z=1
    packing = 'Gamma'
    vasp.kpoints = "Gamma_Mesh\n0\n{}\n{} {} {}".format(packing, x, y, z)
    return vasp

def set_441(vasp: Vasp, structure=None):
    x=4; y=4 ; z=1
    packing = 'Gamma'
    vasp.kpoints = "Gamma_Mesh\n0\n{}\n{} {} {}".format(packing, x, y, z)
    return vasp

def gamma_optimization(vasp: Vasp, structure = None):
    vasp.add_keyword('AUTO_GAMMA', True)
    return vasp

def set_nkred_222(vasp: Vasp, structure=None):
    vasp.add_keyword('nkred', 2)
    return vasp


def set_nkred_221(vasp: Vasp, structure=None):
    vasp.add_keyword('nkredx', 2)
    vasp.add_keyword('nkredy', 2)
    return vasp

def unset_nkred(vasp: Vasp, structure=None):
    vasp.add_keyword('nkredx', 1)
    vasp.add_keyword('nkredy', 1)
    vasp.add_keyword('nkredz', 1)
    vasp.add_keyword('nkred', 1)
    return vasp

######
# TS #
######

def set_dimer(vasp: Vasp, structure=None):
    vasp.add_keyword('ICHAIN', 2)
    return vasp