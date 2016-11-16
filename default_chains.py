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
            assert len(names) == len(workflow)
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

def set_kpoints(x, y, z, packing='Gamma'):
    def kpoints_fxn(vasp, structure=None):
        vasp.kpoints = "Set Mesh\n0\n{}\n{} {} {}".format(packing,x,y,z)
        if packing[0].upper() == 'G':
            vasp.add_keyword('auto_gamma', 'True')
        return vasp
    return kpoints_fxn

def set_spin(magmom):
    def spin_fxn(vasp, structure):
        assert len(magmom) == len(structure)
        vasp.ispin = 2
        vasp.magmom = True
        old_structure = structure.copy
        for _ in range(len(magmom)): # remove all atoms
            structure.pop(0)
        for i in range(len(magmom)): # add atoms again, removing
            x = old_structure[i].pos[0]
            y = old_structure[i].pos[1]
            z = old_structure[i].pos[2]            
            structure.add_atom(x, y, z, structure[i].type, magmom=magmom[i])
        return vasp
    return spin_fxn
