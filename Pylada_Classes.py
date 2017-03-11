#!/usr/bin/env python

from pylada.crystal import read,write
from pymatgen.core import Structure
from Classes_Pymatgen import Poscar
import tempfile

def pylada_to_pmg(pylada_structure):
    '''
    Converts a pylada structure to a pymatgen structure via a temporary file

    :param pylada_structure: Pylada structure to convert
    :return: Structure
    '''
    with tempfile.NamedTemporaryFile() as tmp:
        name = tmp.name
        write.poscar(pylada_structure, name, vasp5=True)
        pymatgen_structure = Poscar.from_file(name).structure # pymatgen.core.Structure
    return pymatgen_structure

def pmg_to_pylada(pmg_structure):
    '''
    Converts a pylada structure to a pymatgen structure via a temporary file

    :param pylada_structure: Pylada structure to convert
    :return: Structure
    '''
    with tempfile.NamedTemporaryFile() as tmp:
        name = tmp.name
        pmg_structure.to('poscar', name)
        pylada_structure = read.poscar(name)
    return pylada_structure
