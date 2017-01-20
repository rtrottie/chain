from default_chains import *
from pylada.vasp.specie import U

def load_optimized_U_species(vasp : Vasp, structure):
    # See vasp/functional.py:  elementName, fileName, max or min oxidation state
    pseudoDir = '$PSEUDO_DIR'
    vasp.add_specie = "Sc", pseudoDir + "/Sc_sv", U("dudarev", "d", 3)    # FERE
    vasp.add_specie = "Ti", pseudoDir + "/Ti_pv", U("dudarev", "d", 4.35) # Wolverton 3+
    vasp.add_specie = "V",  pseudoDir + "/V_pv",  U("dudarev", "d", 4.9)  # Wolverton 2+
    vasp.add_specie = "Cr", pseudoDir + "/Cr_pv", U("dudarev", "d", 3.04) # Wolverton 3+
    vasp.add_specie = "Mn", pseudoDir + "/Mn_pv", U("dudarev", "d", 2.98) # Wolverton 2+
    vasp.add_specie = "Fe", pseudoDir + "/Fe_pv", U("dudarev", "d", 4.04) # Wolverton 2+,
    vasp.add_specie = "Co", pseudoDir + "/Co_pv", U("dudarev", "d", 3.75) # Wolverton 2+
    vasp.add_specie = "Ni",  pseudoDir + "/Ni_pv", U("dudarev", "d", 4.4) # Wolverton 2+
    vasp.add_specie = "Cu", pseudoDir + "/Cu_pv", U("dudarev", "d" , 5.0) # FERE
    vasp.add_specie = "Zn", pseudoDir + "/Zn",                            
    vasp.add_specie = "O",  pseudoDir + "/O"
    vasp.add_specie = "Al", pseudoDir + "/Al"
    vasp.add_specie = "H", pseudoDir + "/H"
    return(vasp)


class WSBulkChain_ferro(CustomChain):
    def __init__(self, vaspobj: Vasp):
        spin = ferro_spin
        standard = [load_default_vasp, ws_standard, ws_bulk, load_optimized_U_species, spin, rough_converge, set_222, set_iopt_7]
        pre_converge   = CustomFunctional(Vasp, standard + [awful_converge, set_gamma, gamma_optimization])
        bad_converge   = CustomFunctional(Vasp, standard + [rough_converge])
        get_nopsin_eig = CustomFunctional(Vasp, standard + [get_eigen_nospin])
        get_eigenvalues= CustomFunctional(Vasp, standard + [get_eigen])
        final_converge = CustomFunctional(Vasp, standard + [full_converge])
        hse            = CustomFunctional(Vasp, standard + [single_point, hse06, set_nkred_222])
        dos            = CustomFunctional(Vasp, standard + [single_point, hse06, set_dos, tetrahedron, all_output])
        names          = ['0_pre_converge', '1_rough_converge', '2_nospin_eig', '3_get_eigenvalues', '4_final_converge', '5_hse']
        return super().__init__([pre_converge, bad_converge, get_nopsin_eig, get_eigenvalues, final_converge, hse], names=names)

class WSBulkChain_anti(CustomChain):
    def __init__(self, vaspobj: Vasp):
        spin = anti_spin
        standard = [load_default_vasp, ws_standard, ws_bulk, load_optimized_U_species, spin, rough_converge, set_222, set_iopt_7]
        pre_converge   = CustomFunctional(Vasp, standard + [awful_converge, set_gamma, gamma_optimization])
        bad_converge   = CustomFunctional(Vasp, standard + [rough_converge])
        get_nopsin_eig = CustomFunctional(Vasp, standard + [get_eigen_nospin, set_algo_normal])
        get_eigenvalues= CustomFunctional(Vasp, standard + [get_eigen])
        final_converge = CustomFunctional(Vasp, standard + [full_converge])
        hse            = CustomFunctional(Vasp, standard + [hse06, set_nkred_222])
        dos            = CustomFunctional(Vasp, standard + [single_point, hse06, tetrahedron, all_output])
        names          = ['0_pre_converge', '1_rough_converge', '2_nospin_eig', '3_get_eigenvalues', '4_final_converge', '5_hse']
        return super().__init__([pre_converge, bad_converge, get_nopsin_eig, get_eigenvalues, final_converge, hse], names=names)

class WSSurfaceChain(CustomChain):
    def __init__(self, vaspobj: Vasp, standard = [], override = []):
        spin = ferro_spin
        standard = [load_default_vasp, ws_standard, ws_surface, load_optimized_U_species, spin, rough_converge, set_221, set_iopt_7, set_kpar_auto] + standard
        pre_converge = CustomFunctional(Vasp, standard + [awful_converge, set_gamma, gamma_optimization, set_algo_fast] + override)
        bad_converge = CustomFunctional(Vasp, standard + [rough_converge, set_algo_fast] + override)
        get_nopsin_eig = CustomFunctional(Vasp, standard + [get_eigen_nospin, set_algo_normal] + override)
        get_eigenvalues = CustomFunctional(Vasp, standard + [get_eigen, set_algo_normal] + override)
        final_converge = CustomFunctional(Vasp, standard + [full_converge, set_algo_fast, all_output ] + override)
        hse = CustomFunctional(Vasp, standard + [single_point, hse06, set_nkred_221, all_output] + override)
#        dos = CustomFunctional(Vasp, standard + [single_point, hse06, set_nkred_221, tetrahedron, all_output, set_dos] + override)
        names          = ['0_pre_converge', '1_rough_converge', '2_nospin_eig', '3_get_eigenvalues', '4_final_converge', '5_hse']
        return super().__init__([pre_converge, bad_converge, get_nopsin_eig, get_eigenvalues, final_converge, hse], names=names, vaspobj=vaspobj)

class WSSurfaceChain_unit(WSSurfaceChain):
    def __init__(self, vaspobj: Vasp):
        return  super().__init__(vaspobj, standard=[set_441])


class WSSurfaceChain_gamma(WSSurfaceChain):
    def __init__(self, vaspobj : Vasp):
        return  super().__init__(vaspobj, standard=[set_gamma, gamma_optimization], override=[unset_nkred])

class WSSurfaceChain_gamma_dimer(WSSurfaceChain):
    def __init__(self, vaspobj : Vasp):
        return  super().__init__(vaspobj, standard=[set_gamma, gamma_optimization, set_dimer], override=[unset_nkred])



spins = {
    'Sc' : 1,
    'Ti' : 2,
    'V'  : 3,
    'Cr' : 4,
    'Mn' : 5,
    'Fe' : 4,
    'Co' : 3,
    'Ni' : 2,
    'Cu' : 1,
    'Zn' : 0,
    'O'  : 0,
    'Al' : 0
}

def ws_standard(vasp: Vasp, structure):
    vasp.isym = 0
    vasp.ismear = -5
    return vasp

def anti_spin(vasp, structure):
    vasp.ispin = 2
    vasp.magmom = '32*0 8*4 8*-4 64*0' 
    return vasp

def ferro_spin(vasp, structure):
    vasp.ispin = 2
    vasp.magmom = True
    for a in structure:
        if a.type in spins:
            a.magmom = spins[a.type]
    return vasp

def ws_bulk(vasp: Vasp, structure):
    vasp.prec = 'Accurate'
    vasp.encut = 500
    vasp.ediff = 1e-6
    vasp.ispin = 2
    vasp.magmom = True
    vasp.magmom = True
    vasp.sigma = 0.01
    vasp.ldau = True
    return vasp

def ws_surface(vasp: Vasp, structure):
    ws_bulk(vasp, structure)
    vasp.encut = 800
    return vasp
