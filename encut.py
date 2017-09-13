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
from pymatgen.core import Structure
from Pylada_Classes import pylada_to_pmg
from math import floor, ceil
import pylada
import os
import math
from Classes_Pymatgen import Poscar

def set_300(vasp: Vasp, structure=None):
    vasp.encut = 300
    return vasp

def set_350(vasp: Vasp, structure=None):
    vasp.encut = 350
    return vasp

def set_400(vasp: Vasp, structure=None):
    vasp.encut = 400
    return vasp

def set_450(vasp: Vasp, structure=None):
    vasp.encut = 450
    return vasp

def set_500(vasp: Vasp, structure=None):
    vasp.encut = 500
    return vasp

def set_550(vasp: Vasp, structure=None):
    vasp.encut = 550
    return vasp

def set_600(vasp: Vasp, structure=None):
    vasp.encut = 600
    return vasp

def set_650(vasp: Vasp, structure=None):
    vasp.encut = 650
    return vasp

def set_700(vasp: Vasp, structure=None):
    vasp.encut = 700
    return vasp

def set_750(vasp: Vasp, structure=None):
    vasp.encut = 750
    return vasp

def set_800(vasp: Vasp, structure=None):
    vasp.encut = 800
    return vasp

def set_850(vasp: Vasp, structure=None):
    vasp.encut = 850
    return vasp

def set_900(vasp: Vasp, structure=None):
    vasp.encut = 900
    return vasp