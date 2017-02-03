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
from default_chains import CustomChain
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

    def get_bandgap_from_aexx(self, aexx):
        vasprun_location = os.path.join(str(aexx).zfill(2), self.names[-1], 'vasprun.xml')
        vasprun = Vasprun(vasprun_location, parse_dos=False, parse_eigen=False, parse_projected_eigen=False, parse_potcar_file=False)
        vasprun.get_band_structure().get_band_gap()

    def find_aexx(self, aexx_low : int, aexx_high : int):
        '''
        Does a binary search to find correct value of AEXX
        :param aexx_low:
        :param aexx_high:
        :return:
        '''
        aexx_low = max(0, aexx_low)
        aexx_high = min(99, aexx_high)

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