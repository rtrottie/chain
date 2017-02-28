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
from pylada.crystal import read,write
from pymatgen.io.vasp import Vasprun
from math import floor
import pylada
import os
import math
from Classes_Pymatgen import Poscar


class CustomFunctional(object):
    def __init__(self, base: Vasp, modifications: list, type='relax'):
        self.base = base
        self.modifications = modifications
        self.type = type
        return

class CustomChain(object):
    def __init__(self, functionals: list, names=None, vaspobj:Vasp=None, basename='', initial_structure=None, final_structure=None):
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
        self.initial_structure = initial_structure
        self.final_strucutre = final_structure
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
        '''
        The important function run when __call__ is invoked
        :param name: name (subfolder) of current run
        :param workflow: The workflow to be run
        :param structure: Structure that will be relaxed
        :param outdir:  Directory for output
        :param previous: Previous directory (for initialization)
        :param kwargs:
        :return:
        '''
        vasp = workflow.base(copy=deepcopy(self.vasp))
        structure_ = structure.copy()
        outdir = os.getcwd() if outdir is None else RelativePath(outdir).path
        for modification in workflow.modifications:
            vasp = modification(vasp, structure_)
        ## if this calculation has not been done run it
        params = deepcopy(kwargs)
        fulldir = os.path.join(outdir, name)
        if workflow.type == 'neb':  #TODO:  Support multiple images
            images = vasp.images+2 # for initial and final point
            for image in range(images):
                image_dir = os.path.join(fulldir, str(image).zfill(2))
                os.makedirs(image_dir, exist_ok=True)
                if image == 0:
                    write.poscar(self.initial_structure, os.path.join(image_dir, 'POSCAR'))
                elif image == len(images)-1:
                    write.poscar(self.final_structure, os.path.join(image_dir, 'POSCAR'))
                elif image == floor(len(images)/2):
                    write.poscar(self.structure_, os.path.join(image_dir, 'POSCAR'))


        output = vasp(structure_, outdir=fulldir, restart=previous ,**params)

        if not output.success:
            raise ExternalRunFailed("VASP calculation did not succeed.")
        return output

    # Creating the workflow
    def call_with_output(self, structure, outdir=None, names=None, functionals=None, previous=None, **kwargs ):
        if not names:
            names = self.names
        if not functionals:
            functionals = self.functionals
        # make this function stateless.
        for i in range(len(functionals)):
            name = names[i]
            workflow = functionals[i]
            print(previous)
            previous = self.run_calculation(name, workflow, structure, outdir, previous, **kwargs)
        fulldir = os.path.join(outdir, name)
        return (self.Extract(fulldir), previous)


    def __call__(self, structure, outdir=None, names=None, functionals=None, previous=None, **kwargs ):
        (extract, _) = self.call_with_output(structure, outdir=outdir, names=names, functionals=functionals, previous=previous)
        return extract

class OptimizedParametersChain(CustomChain):
    def __init__(self, functionals: list, bandgap:float=None, names=None, vaspobj:Vasp=None, basename=''):
        '''
        Runs a series of workflows
        Args:
            functionals (list): list of lists.  Each sublist should be a series of functions that take, and statically modify, a pylada.vasp.relax.Relax object.
            nupdown_functionals (list) : list of lists.  Each sublist should be a series of functions that take, and statically modify, a pylada.vasp.relax.Relax object.
            nupdowns (int list): nupdowns that will be considered
            vaspobj: TODO
        '''
        self.bandgap =  bandgap
        return  super().__init__(functionals, names=names, vaspobj=vaspobj, basename=basename)


    def get_bandgap_from_aexx(self, structure, aexx, outdir=None, previous=None):
        vasprun_location = os.path.join(outdir, str(aexx).zfill(2), self.names[-1], 'vasprun.xml')
        # try:
        #     vasprun = Vasprun(vasprun_location, parse_projected_eigen=False)
        #     band_gap = vasprun.get_band_structure().get_band_gap()['energy']
        # except:
        def set_aexx(vasp: Vasp, structure=None):
            vasp.add_keyword('aexx', aexx/100)
            return vasp
        for x in self.functionals: # Set nupdown
            x.modifications.append(set_aexx)
        (_, output) = super().call_with_output(structure, outdir=os.path.join(outdir, str(aexx).zfill(2)), previous=previous)
        vasprun = Vasprun(vasprun_location, parse_projected_eigen=False)
        band_gap = vasprun.get_band_structure().get_band_gap()['energy']
        return (band_gap, output)
    def find_aexx(self, structure, aexx_low : int, aexx_high : int, outdir=None, previous=None):
        '''
        Does a binary search to find correct value of AEXX
        :param aexx_low:
        :param aexx_high:
        :return:
        '''
        aexx_center = int((aexx_high + aexx_low) / 2)

        (bg_low, _)  = self.get_bandgap_from_aexx(structure, aexx_low, outdir, previous=previous)
        (bg_high, _) = self.get_bandgap_from_aexx(structure, aexx_high, outdir, previous=previous)
        (bg_center, previous) = self.get_bandgap_from_aexx(structure, aexx_center, outdir, previous=previous)
        if aexx_high - aexx_low <= 1:
            if abs(bg_low - self.bandgap) <= abs(bg_high - self.bandgap):
                return aexx_low
            else:
                return aexx_high
        elif bg_low > self.bandgap:
            raise Exception('Bandgap low above desired bandgap')
        elif bg_high < self.bandgap:
            raise Exception('Bandgap high below desired bandgap')
        elif bg_center > self.bandgap:
            return self.find_aexx(structure, aexx_low, aexx_center, outdir, previous=previous)
        elif bg_center < self.bandgap:
            return self.find_aexx(structure, aexx_center, aexx_high, outdir, previous=previous)

    def get_energy_from_encut(self, structure, encut, outdir=None, convergence_value=1e-4, previous=None):
        vasprun_location = os.path.join(outdir, str(encut).zfill(4), self.names[-1], 'vasprun.xml')
        # try:
        #     vasprun = Vasprun(vasprun_location, parse_projected_eigen=False)
        #     energy = vasprun.final_energy
        # except:
        def set_encut(vasp: Vasp, structure=None):
            vasp.encut = encut
            vasp.ediff = convergence_value/1000
            return vasp
        for x in self.functionals: # Set nupdown
            x.modifications.append(set_encut)
        (_,output) = super().call_with_output(structure, outdir=os.path.join(outdir, str(encut).zfill(4)), previous=previous)
        vasprun = Vasprun(vasprun_location, parse_projected_eigen=False)
        energy = vasprun.final_energy / vasprun.final_structure.num_sites
        return (energy, output)
    def get_encut(self, structure, encut_low : int, encut_high : int, optimal_energy : float, convergence_value : float,  outdir : str, previous=None):
        '''
        Does a binary search to find optimal encut value.  Converges to value within convergence_value variable
        Starts off incrementing encut by 500 to find asymptote

        :param structure: Structure to find optimal encut value for
        :param encut_low: Lower bound for ENCUT
        :param encut_high: Upper bound for ENCUT
        :param optimal_energy: Determined optimal value for energy
        :param convergence_value: value to converge within (eV/atom)
        :param outdir: outdir to write to
        :return: :type int
        '''
        assymptote_increment = 250
        encut_increment = 25
        def encut_round(i: int):
            return round(i/encut_increment)*encut_increment
        (energy_low, _)  = self.get_energy_from_encut(structure, encut_low , outdir, convergence_value, previous=previous)
        (energy_high, output) = self.get_energy_from_encut(structure, encut_high, outdir, convergence_value, previous=previous)

        if optimal_energy >= 0:  # haven't reached assymptote
            if abs(energy_high - energy_low) <= convergence_value:  # reached asymptote
                return self.get_encut(structure, encut_low-assymptote_increment, encut_low+encut_increment, energy_high, convergence_value, outdir, previous=output)
            else:  # keep searching
                return self.get_encut(structure, encut_high, encut_high + assymptote_increment, 0, convergence_value, outdir, previous=output)
        else:  # Do Binary search
            encut_center = encut_round((encut_high+encut_low)/2)
            (energy_center,_) = self.get_energy_from_encut(structure, encut_center, outdir, convergence_value, previous=output)
            if (encut_high-encut_low) == encut_increment: # Binary search is at center
                return encut_high
            elif abs(energy_center - optimal_energy) <= convergence_value:  # reached convergence
                return self.get_encut(structure, encut_low, encut_center, optimal_energy, convergence_value, outdir, previous=output)
            else: # center not converged
                return self.get_encut(structure, encut_center, encut_high, optimal_energy, convergence_value, outdir, previous=output)

    def get_energy_from_kpoint(self, structure, kpoint, outdir=None, convergence_value=1e-4, previous=None):
        folder = '{0}'.format(str(kpoint).zfill(2))
        vasprun_location = os.path.join(outdir, folder, self.names[-1], 'vasprun.xml')
        # try:
        #     vasprun = Vasprun(vasprun_location, parse_projected_eigen=False)
        #     energy = vasprun.final_energy
        # except:
        def set_kpoint(vasp: Vasp, structure):
            lengths = [sum([x ** 2 for x in structure.cell.transpose()[i]]) ** (1 / 2) for i in range(3)] # using distance formula to get vector lengths
            kpoints = [math.ceil(min(lengths) / x * kpoint) for x in lengths] # scaling number of kpoints for shorter vectors

            packing = 'Gamma'
            vasp.kpoints = "Gamma_Mesh\n0\n{}\n{} {} {}".format(packing, kpoints[0], kpoints[1], kpoints[2])
            vasp.ediff = convergence_value/1000
            return vasp
        for x in self.functionals: # Set nupdown
            x.modifications.append(set_kpoint)
        (_,output) = super().call_with_output(structure, outdir=os.path.join(outdir, folder), previous=previous)
        vasprun = Vasprun(vasprun_location, parse_projected_eigen=False)
        energy = vasprun.final_energy / vasprun.final_structure.num_sites
        return (energy, output)
    def get_kpoints(self, structure, kpoint, convergence_value: float, outdir: str, previous=None):
        (energy_asymptote, output) = self.get_energy_from_kpoint(structure, kpoint, outdir, convergence_value, previous=previous)
        if kpoint <= 5:
            return self.get_kpoints(structure, kpoint+1, convergence_value, outdir, previous=output)
        else:
            (energy_low,_) = self.get_energy_from_kpoint(structure, kpoint-3, outdir, convergence_value, previous=output)
            (energy_asymptote_low,_) = self.get_energy_from_kpoint(structure, kpoint-1, outdir, convergence_value, previous=output)
            if abs(energy_low - energy_asymptote) <= convergence_value and abs(energy_low - energy_asymptote_low) <= convergence_value:
                return kpoint - 3
            else:
                return self.get_kpoints(structure, kpoint + 1, convergence_value, outdir, previous=output)



    def __call__(self, structure, outdir=None, **kwargs):
        kpoint = self.get_kpoints(structure, 3, 0.0005, outdir=os.path.join(outdir, 'get_kpoints'))
        def set_kpoint(vasp: Vasp, structure=None):
            packing = 'Gamma'
            vasp.kpoints = "Gamma_Mesh\n0\n{0}\n{1} {1} {1}".format(packing, kpoint)
            return vasp
        for x in self.functionals: # Set nupdown
            x.modifications.append(set_kpoint)

        encut = self.get_encut(structure, 300, 800, 0, 0.0005 ,outdir=os.path.join(outdir, 'get_encut'))
        def set_encut(vasp: Vasp, structure=None):
            vasp.encut = encut
            return vasp
        for x in self.functionals: # Set nupdown
            x.modifications.append(set_encut)

        aexx = self.find_aexx(structure, 0, 99, outdir=os.path.join(outdir, 'get_aexx'))  / 100
        def set_aexx(vasp: Vasp, structure=None):
            vasp.encut = encut
            return vasp
        for x in self.functionals: # Set nupdown
            x.modifications.append(set_aexx)

        with open(os.path.join(outdir, 'INCAR.defaults'), 'w') as f:
            f.write('AEXX = {}\n'.format(aexx))
            f.write('ENCUT = {}\n'.format(encut))
            f.write('KPOINTS = {}\n'.format(kpoint))
        return super().__call__(structure, outdir=outdir)

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
        if energies:
            nup = min(energies.keys(), key=(lambda key: energies[key]) )
            def set_nupdown(vasp: Vasp, structure=None):
                vasp.nupdown = nup
                return vasp
            for x in self.functionals: # Set nupdown
                x.modifications.append(set_nupdown)
        return super().__call__(structure, outdir=outdir, **kwargs)

class TSSpinCustomChain(SpinCustomChain):
    def __init__(self, functionals: list, initial_structure: Structure, final_structure: Structure = None, names=None, vaspobj:Vasp=None, basename=''):
        self.initial = initial_structure
        self.final = final_structure
        return super().__init__(functionals, names=names, vaspobj=vaspobj, basename=basename)

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
    vasp.add_keyword('kpar', 1)
    vasp.npar = int(nodes)
    return vasp

def set_kpar_per_node(vasp: Vasp, structure=None):
    nodes = int(os.environ['PBS_NUM_NODES'])
    vasp.add_keyword('kpar', nodes)
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

def no_hse06(vasp: Vasp, structure=None):
    vasp.add_keyword('lhfcalc', False)
    vasp.ismear =  -5
    vasp.algo = 'Normal'
    vasp.ldau = True
    return vasp

def set_dos(vasp: Vasp, structure=None):
    vasp.algo = 'None'
    return vasp

def tetrahedron(vasp: Vasp, structure=None):
    vasp.ismear = -5
    return vasp

def set_npar_2(vasp: Vasp, structure=None):
    vasp.npar = 2
    return vasp

##########
# IONIC  #
##########

def set_iopt_7(vasp: Vasp, structure=None):
    vasp.ibrion = 3
    vasp.add_keyword('potim',  0)
    vasp.add_keyword('iopt', 7)
    return vasp

def single_point(vasp: Vasp, structure=None):
    vasp.istart = 1
    vasp.icharg = 1
    vasp.ibrion = -1
    vasp.nelmin = 3
    vasp.nsw = 0
    vasp.ediff = 1e-5
    vasp.add_keyword('lmaxmix', None)
    vasp.add_keyword('iopt', 0)
    return vasp

def cell_relax(vasp: Vasp, structure=None):
    vasp.isif = 3
    vasp.ibrion = 1
    vasp.potim = 0.4
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

def set_prec_normal(vasp: Vasp, structure=None):
    vasp.prec = 'normal'
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

def set_333(vasp: Vasp, structure=None):
    x=3; y=3 ; z=3
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
    vasp.add_keyword('auto_gamma', True)
    return vasp

def set_nkred_222(vasp: Vasp, structure=None):
    vasp.add_keyword('nkred', 2)
    return vasp

def set_nkred_333(vasp: Vasp, structure=None):
    vasp.add_keyword('nkred', 3)
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
    vasp.add_keyword('ichain', 2)
    vasp.add_keyword('images', 0)
    return vasp

def set_single_neb(vasp: Vasp, structure=None):
    vasp.add_keyword('ichain', 0)
    vasp.add_keyword('images', 1)
    vasp.add_keyword('lclimb', True)
    return vasp
