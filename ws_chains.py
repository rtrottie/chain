from default_chains import *
from pylada.vasp.specie import U

def load_optimized_U_species(vasp, structure):
    # See vasp/functional.py:  elementName, fileName, max or min oxidation state
    pseudoDir = '$PSEUDO_DIR'
    vasp.add_specie = "Sc", pseudoDir + "/Sc_sv", U("dudarev", "d", 3)    # FERE
    vasp.add_specie = "Ti", pseudoDir + "/Ti_pv", U("dudarev", "d", 4.35) # Wolverton 3+
    vasp.add_specie = "V",  pseudoDir + "/V_pv",  U("dudarev", "d", 4.9)  # Wolverton 2+
    vasp.add_specie = "Cr", pseudoDir + "/Cr_pv", U("dudarev", "d", 3.04) # Wolverton 3+
    vasp.add_specie = "Mn", pseudoDir + "/Mn_pv", U("dudarev", "d", 2.98) # Wolverton 2+
    vasp.add_specie = "Fe", pseudoDir + "/Fe_pv", U("dudarev", "d", 4.04) # Wolverton 2+
    vasp.add_specie = "Co", pseudoDir + "/Co_pv", U("dudarev", "d", 3.75) # Wolverton 2+
    vasp.add_specie = "Ni",  pseudoDir + "/V_pv", U("dudarev", "d", 4.4)  # Wolverton 2+
    vasp.add_specie = "Cu", pseudoDir + "/Cu_pv", U("dudarev", "d" , 5.0) # FERE
    vasp.add_specie = "Zn", pseudoDir + "/Zn",                            
    vasp.add_specie = "O",  pseudoDir + "/O"
    vasp.add_specie = "Al", pseudoDir + "/Al"
    vasp.add_specie = "H", pseudoDir + "/H"
    return(vasp)

class WSBulkChain_ferro(CustomChain):
    def __init__(self, vaspobj):
        spin = ferro_spin
        pre_converge = [load_default_vasp, load_optimized_U_species, spin, awful_converge, set_gamma, set_iopt_7]
        bad_converge = [load_default_vasp, load_optimized_U_species, spin, rough_converge, set_222, set_iopt_7]
        get_eigenvalues = [load_default_vasp, load_optimized_U_species, spin, get_eigen, set_222, set_iopt_7]
        final_converge = [load_default_vasp, load_optimized_U_species, spin, full_converge, set_222, set_iopt_7]
        hse = [load_default_vasp, load_optimized_U_species, spin, single_point, hse06, set_222, set_iopt_7]
        dos = [load_default_vasp, load_optimized_U_species, spin, single_point, hse06, set_dos, set_222, tetrahedron, set_iopt_7]
        names = ['0_pre_converge', '1_rough_converge', '2_get_eigenvalues', '3_final_converge', '4_hse', '5_dos']
        return super().__init__([pre_converge, bad_converge, get_eigenvalues, final_converge, hse, dos], names=names, vaspobj=vaspobj)

def anti_spin(vasp, structure):
    vasp.ispin = 2
    vasp.magmom = '32*0 8*4 8*-4 64*0' 
    return vasp

def ferro_spin(vasp, structure):
    vasp.ispin = 2
    vasp.magmom = '32*0 16*4 64*0' 
    return vasp

