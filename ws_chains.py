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
    vasp.add_specie = "Fe", pseudoDir + "/Fe_pv", U("dudarev", "d", 4.04) # Wolverton 2+,
    vasp.species['Fe'].magmom = 4
    vasp.add_specie = "Co", pseudoDir + "/Co_pv", U("dudarev", "d", 3.75) # Wolverton 2+
    vasp.add_specie = "Ni",  pseudoDir + "/V_pv", U("dudarev", "d", 4.4)  # Wolverton 2+
    vasp.add_specie = "Cu", pseudoDir + "/Cu_pv", U("dudarev", "d" , 5.0) # FERE
    vasp.add_specie = "Zn", pseudoDir + "/Zn",                            
    vasp.add_specie = "O",  pseudoDir + "/O"
    vasp.add_specie = "Al", pseudoDir + "/Al"
    vasp.add_specie = "H", pseudoDir + "/H"
    return(vasp)


class WSBulkChain_ferro(CustomChain):
    def __init__(self, vaspobj: Vasp):
        spin = ferro_spin
        standard = [ws_bulk, load_default_vasp, load_optimized_U_species, spin, rough_converge, set_222, set_iopt_7]
        pre_converge   = CustomFunctional(Vasp, standard + [awful_converge, set_gamma, gamma_optimization])
        bad_converge   = CustomFunctional(Vasp, standard + [rough_converge])
        get_nopsin_eig = CustomFunctional(Vasp, standard + [get_eigen_nospin])
        get_eigenvalues= CustomFunctional(Vasp, standard + [get_eigen])
        final_converge = CustomFunctional(Vasp, standard + [full_converge])
        hse            = CustomFunctional(Vasp, standard + [single_point, hse06])
        dos            = CustomFunctional(Vasp, standard + [single_point, hse06, set_dos, tetrahedron])
        names          = ['0_pre_converge', '1_rough_converge', '2_nospin_eig', '3_get_eigenvalues', '4_final_converge', '5_hse', '6_dos']
        return super().__init__([pre_converge, bad_converge, get_nopsin_eig, get_eigenvalues, final_converge, hse, dos], names=names)

class WSBulkChain_anti(CustomChain):
    def __init__(self, vaspobj: Vasp):
        spin = anti_spin
        standard = [ws_bulk, load_default_vasp, load_optimized_U_species, spin, rough_converge, set_222, set_iopt_7]
        pre_converge   = CustomFunctional(Vasp, standard + [awful_converge, set_gamma, gamma_optimization])
        bad_converge   = CustomFunctional(Vasp, standard + [rough_converge])
        get_eigenvalues= CustomFunctional(Vasp, standard + [get_eigen])
        final_converge = CustomFunctional(Vasp, standard + [full_converge])
        hse            = CustomFunctional(Vasp, standard + [hse06])
        dos            = CustomFunctional(Vasp, standard + [single_point, hse06, set_dos, tetrahedron])
        names          = ['0_pre_converge', '1_rough_converge', '2_get_eigenvalues', '3_final_converge', '4_hse', '5_dos']
        return super().__init__([pre_converge, bad_converge, get_eigenvalues, final_converge, hse, dos], names=names)

class WSSurfaceChain(CustomChain):
    def __init__(self, vaspobj: Vasp, standard = []):
        spin = ferro_spin
        standard = [ws_surface, load_default_vasp, load_optimized_U_species, rough_converge, set_221, set_iopt_7] + standard
        pre_converge = CustomFunctional(Vasp, standard + [awful_converge, set_gamma, gamma_optimization, set_algo_fast])
        bad_converge = CustomFunctional(Vasp, standard + [rough_converge, set_algo_fast])
        get_nopsin_eig = CustomFunctional(Vasp, standard + [get_eigen_nospin, set_algo_normal])
        get_eigenvalues = CustomFunctional(Vasp, standard + [get_eigen, set_algo_normal])
        final_converge = CustomFunctional(Vasp, standard + [full_converge, set_algo_fast])
        hse = CustomFunctional(Vasp, standard + [single_point, hse06])
        dos = CustomFunctional(Vasp, standard + [single_point, hse06, set_dos, tetrahedron])
        names          = ['0_pre_converge', '1_rough_converge', '2_nospin_eig', '3_get_eigenvalues', '4_final_converge', '5_hse', '6_dos']
        return super().__init__([pre_converge, bad_converge, get_nopsin_eig, get_eigenvalues, final_converge, hse, dos], names=names)

class WSSurfaceChain_unit(WSSurfaceChain):
    def __init__(self, vaspobj: Vasp):
        return  super().__init__(vaspobj, standard=[set_441])

def anti_spin(vasp, structure):
    vasp.ispin = 2
    vasp.magmom = '32*0 8*4 8*-4 64*0' 
    return vasp

def ferro_spin(vasp, structure):
    vasp.ispin = 2
    vasp.magmom = '32*0 16*4 64*0' 
    return vasp

def ws_bulk(vasp: Vasp, structure):
    vasp.prec = 'Accurate'
    vasp.encut = 500
    vasp.ediff = 1e-6
    vasp.ispin = 2
    vasp.ismear = 0
    vasp.sigma = 0.01
    vasp.ldau = True
    return vasp

def ws_surface(vasp: Vasp, structure):
    ws_bulk(vasp, structure)
    vasp.encut = 800
    return vasp