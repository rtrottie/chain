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

    def get_bandgap_from_aexx(self, structure, aexx):
        vasprun_location = os.path.join(str(aexx).zfill(2), self.names[-1], 'vasprun.xml')
        try:
            vasprun = Vasprun(vasprun_location, parse_projected_eigen=False)
            band_gap = vasprun.get_band_structure().get_band_gap()['energy']
        except:
            def set_aexx(vasp: Vasp, structure=None):
                vasp.add_keyword('AEXX', aexx/100)
                return vasp
            for x in self.functionals: # Set nupdown
                x.modifications.append(set_aexx)
            super().__call__(structure=structure, outdir=str(aexx).zfill(2))
            vasprun = Vasprun(vasprun_location, parse_projected_eigen=False)
            band_gap = vasprun.get_band_structure().get_band_gap()['energy']
        return band_gap

    def find_aexx(self, structure, aexx_low : int, aexx_high : int):
        '''
        Does a binary search to find correct value of AEXX
        :param aexx_low:
        :param aexx_high:
        :return:
        '''
        aexx_center = int((aexx_high + aexx_low) / 2)

        bg_low  = self.get_bandgap_from_aexx(structure, aexx_low)
        bg_high = self.get_bandgap_from_aexx(structure, aexx_high)
        bg_center = self.get_bandgap_from_aexx(structure, aexx_center)
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
            self.find_aexx(structure, aexx_low, aexx_center)
        elif bg_center < self.bandgap:
            self.find_aexx(structure, aexx_center, aexx_high)


    def __call__(self, structure, outdir=None, **kwargs):
        return self.find_aexx(structure, 0, 99)

class BulkHSE(AEXX):
    def __init__(self, vaspobj: Vasp(), nupdowns, standard=[], override=[], final_step='5_hse' ):
        standard = [load_default_vasp, cell_relax, set_iopt_7, bulk_standard]
        pbe = CustomFunctional(vaspobj, standard)
        hse = CustomFunctional(Vasp, standard + [hse06])
        hse_single = CustomFunctional(Vasp, standard + [hse06, single_point, all_output])
        names = ['0_pbe', '1_hse', '2_hse_singlepoint']
        super().__init__([pbe, hse, hse_single], names=names, vaspobj=vaspobj)

def bulk_standard(vasp: Vasp, structure):
    # Start
    vasp.istart = 0
    vasp.icharg = 2
    # Electronic
    vasp.add_keyword('GGA', 'PS')
    vasp.isym = 0
    vasp.ismear = -5
    vasp.prec = "Accurate"
    vasp.nelm = 60
    vasp.ediff = 1e-5
    vasp.nelmdl = 0
    # Ionic
    vasp.nsw = 5000
    vasp.ediffg = -0.005
    # Output
    vasp.lwave = True
    vasp.lcharg = True

    # TODO: Get these automatically
    x=4; y=4 ; z=4
    packing = 'Gamma'
    vasp.kpoints = "Gamma_Mesh\n0\n{}\n{} {} {}".format(packing, x, y, z)
    vasp.encut = 500
    return vasp