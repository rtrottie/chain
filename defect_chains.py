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
from default_chains import *
from pymatgen.io.vasp.outputs import Vasprun
from Classes_Pymatgen import Incar

class AEXX(CustomChain):
    def __init__(self, functionals: list, bandgap : float, names=None, vaspobj:Vasp=None, basename=''):
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

    def get_bandgap_from_aexx(self, structure, aexx, outdir=None):
        vasprun_location = os.path.join(outdir, str(aexx).zfill(2), self.names[-1], 'vasprun.xml')
        try:
            vasprun = Vasprun(vasprun_location, parse_projected_eigen=False)
            band_gap = vasprun.get_band_structure().get_band_gap()['energy']
        except:
            def set_aexx(vasp: Vasp, structure=None):
                vasp.add_keyword('AEXX', aexx/100)
                return vasp
            for x in self.functionals: # Set nupdown
                x.modifications.append(set_aexx)
            super().__call__(structure, outdir=os.path.join(outdir, str(aexx).zfill(2)))
            vasprun = Vasprun(vasprun_location, parse_projected_eigen=False)
            band_gap = vasprun.get_band_structure().get_band_gap()['energy']
        return band_gap

    def find_aexx(self, structure, aexx_low : int, aexx_high : int, outdir=None):
        '''
        Does a binary search to find correct value of AEXX
        :param aexx_low:
        :param aexx_high:
        :return:
        '''
        aexx_center = int((aexx_high + aexx_low) / 2)

        bg_low  = self.get_bandgap_from_aexx(structure, aexx_low, outdir)
        bg_high = self.get_bandgap_from_aexx(structure, aexx_high, outdir)
        bg_center = self.get_bandgap_from_aexx(structure, aexx_center, outdir)
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
            self.find_aexx(structure, aexx_low, aexx_center, outdir)
        elif bg_center < self.bandgap:
            self.find_aexx(structure, aexx_center, aexx_high, outdir)


    def __call__(self, structure, outdir=None, **kwargs):
        aexx = self.find_aexx(structure, 0, 99, outdir=os.path.join(outdir, 'get_aexx'))
        with open(os.path.join(outdir, 'INCAR.aexx'), 'w') as f:
            f.write('AEXX = {}'.format(aexx))
        return aexx


class BulkPBE(OptimizedParametersChain):
    def __init__(self, vaspobj: Vasp, bandgap:float=None, standard=[], override=[], final_step='5_hse' ):
        standard = [load_default_vasp, cell_relax, bulk_standard, load_species, set_kpar_per_node, set_npar_2]
        pbe = CustomFunctional(Vasp, standard)
        pbe_single = CustomFunctional(Vasp, standard + [single_point, all_output])
        names = ['1_pbe', '2_pbe_singlepoint']
        super().__init__([pbe, pbe_single], bandgap=bandgap, names=names, vaspobj=vaspobj)

class BulkHSE(OptimizedParametersChain):
    def __init__(self, vaspobj: Vasp, bandgap:float=None, standard=[], override=[], final_step='5_hse' ):
        standard = [load_default_vasp, cell_relax, bulk_standard, load_species, set_npar_2, set_kpar_per_node]
        # pbe = CustomFunctional(Vasp, standard)
        hse = CustomFunctional(Vasp, standard + [hse06, set_npar_2])
        hse_single = CustomFunctional(Vasp, standard + [hse06, single_point, all_output, set_npar_2])
        names = ['1_hse', '2_hse_singlepoint']
        super().__init__([hse, hse_single], bandgap=bandgap, names=names, vaspobj=vaspobj)

class DefectHSE(CustomChain):
    def __init__(self, vaspobj: Vasp, standard=[], override=[], final_step='5_hse' ):
        standard = [load_default_vasp, bulk_standard, load_species, set_333, set_nkred_333, set_kpar_auto, set_iopt_7, set_isym_0] + standard
        pbe = CustomFunctional(Vasp, standard)
        pbe_single = CustomFunctional(Vasp, standard + [single_point, all_output])
        hse = CustomFunctional(Vasp, standard + [hse06])
        hse_single = CustomFunctional(Vasp, standard + [hse06, single_point, all_output])
        names = ['1_pbe', '2_pbe_singlepoint', '3_hse', '4_hse_singlepoint']
        super().__init__([pbe, pbe_single, hse, hse_single], names=names, vaspobj=vaspobj )

class DefectHSE_large(CustomChain):
    def __init__(self, vaspobj: Vasp, standard=[], override=[], final_step='5_hse' ):
        standard = [load_default_vasp, bulk_standard, load_species, set_222, set_nkred_222, set_kpar_auto, set_iopt_7, set_isym_0] + standard
        pbe = CustomFunctional(Vasp, standard)
        pbe_single = CustomFunctional(Vasp, standard + [single_point, all_output])
        hse = CustomFunctional(Vasp, standard + [hse06])
        hse_single = CustomFunctional(Vasp, standard + [hse06, single_point, all_output])
        names = ['1_pbe', '2_pbe_singlepoint', '3_hse', '4_hse_singlepoint']
        super().__init__([pbe, pbe_single, hse, hse_single], names=names, vaspobj=vaspobj )

def bulk_standard(vasp: Vasp, structure):
    # Electronic
    vasp.add_keyword('GGA', 'PS')
    vasp.ismear = -5
    vasp.prec = "Accurate"
    vasp.nelm = 60
    vasp.ediff = 1e-5
    vasp.nelmdl = 0
    # Ionic
    vasp.nsw = 5000
    vasp.ediffg = -0.02
    # Output
    vasp.lwave = True
    vasp.lcharg = True

    #TODO: Get these automatically
    x=2; y=2 ; z=2
    packing = 'Gamma'
    vasp.kpoints = "Gamma_Mesh\n0\n{}\n{} {} {}".format(packing, x, y, z)
    return vasp

def load_species(vasp : Vasp, structure):
    # See vasp/functional.py:  elementName, fileName, max or min oxidation state
    pseudoDir = '$PSEUDO_DIR'
    vasp.add_specie = "Zn", pseudoDir + "/Zn"
    vasp.add_specie = "O",  pseudoDir + "/O"
    vasp.add_specie = "Al", pseudoDir + "/Al"
    vasp.add_specie = "Cd", pseudoDir + "/Cd"
    vasp.add_specie = "Te", pseudoDir + "/Te"
    vasp.add_specie = "Ga", pseudoDir + "/Ga"
    vasp.add_specie = "As", pseudoDir + "/As"
    vasp.add_specie = "In", pseudoDir + "/In"
    vasp.add_specie = "P", pseudoDir + "/P"
    vasp.add_specie = "Bi", pseudoDir + "/Bi_d"
    return(vasp)