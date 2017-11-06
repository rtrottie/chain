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
    
    vasp.add_specie = "Sr", pseudoDir + "/Sr_sv"
    vasp.add_specie = "Ba", pseudoDir + "/Ba_sv"
    vasp.add_specie = "La", pseudoDir + "/La" # TODO Determine U
    return(vasp)

class WSBulkChain(SpinCustomChain):
    def __init__(self, vaspobj: Vasp(), nupdowns, standard=[], override=[], final_step='5_hse' ):
        standard = [load_default_vasp, ws_standard, ws_bulk, load_optimized_U_species, rough_converge, set_222, set_iopt_7]
        gamma = [set_gamma, gamma_optimization]
        pre_converge   = CustomFunctional(Vasp, standard + [awful_converge, set_algo_fast] + gamma + override)
        bad_converge   = CustomFunctional(Vasp, standard + [rough_converge, set_algo_fast] + override)
        bad_converge_gamma   = CustomFunctional(Vasp, standard + [rough_converge, set_algo_fast] + gamma + override)
        get_nopsin_eig = CustomFunctional(Vasp, standard + [get_eigen_nospin] + override)
        get_nopsin_eig_gamma = CustomFunctional(Vasp, standard + [get_eigen_nospin] + gamma + override)
        get_eigenvalues= CustomFunctional(Vasp, standard + [get_eigen] + override)
        get_eigenvalues_gamma = CustomFunctional(Vasp, standard + [get_eigen] + gamma + override)
        final_converge = CustomFunctional(Vasp, standard + [full_converge] + override)
        final_converge_gamma = CustomFunctional(Vasp, standard + [full_converge] + gamma + override)
        hse            = CustomFunctional(Vasp, standard + [hse06, set_algo_damp, set_nkred_222] + override + [single_point])
        names          = ['0_pre_converge', '1_rough_converge', '2_nospin_eig', '3_get_eigenvalues', '4_final_converge', final_step]
        nupdown_functionals = [pre_converge, bad_converge_gamma, get_nopsin_eig_gamma, get_eigenvalues_gamma, final_converge_gamma]
        super().__init__([pre_converge, bad_converge, get_nopsin_eig, get_eigenvalues, final_converge, hse],
                         nupdown_functionals=nupdown_functionals, nupdowns=nupdowns, names=names, vaspobj=vaspobj)

class WSBulkChain_small(SpinCustomChain):
    def __init__(self, vaspobj: Vasp(), nupdowns, standard=[], override=[], final_step='5_hse' ):
        standard = [load_default_vasp, ws_standard, ws_bulk, load_optimized_U_species, rough_converge, set_444, set_iopt_7]
        gamma = [set_gamma, gamma_optimization]
        pre_converge   = CustomFunctional(Vasp, standard + [awful_converge, set_algo_fast] + gamma + override)
        bad_converge   = CustomFunctional(Vasp, standard + [rough_converge, set_algo_fast] + override)
        bad_converge_gamma   = CustomFunctional(Vasp, standard + [rough_converge, set_algo_fast] + gamma + override)
        get_nopsin_eig = CustomFunctional(Vasp, standard + [get_eigen_nospin] + override)
        get_nopsin_eig_gamma = CustomFunctional(Vasp, standard + [get_eigen_nospin] + gamma + override)
        get_eigenvalues= CustomFunctional(Vasp, standard + [get_eigen] + override)
        get_eigenvalues_gamma = CustomFunctional(Vasp, standard + [get_eigen] + gamma + override)
        final_converge = CustomFunctional(Vasp, standard + [full_converge] + override)
        final_converge_gamma = CustomFunctional(Vasp, standard + [full_converge] + gamma + override)
        names          = ['0_pre_converge', '1_rough_converge', '2_nospin_eig', '3_get_eigenvalues', '4_final_converge']
        nupdown_functionals = [pre_converge, bad_converge_gamma, get_nopsin_eig_gamma, get_eigenvalues_gamma, final_converge_gamma]
        super().__init__([pre_converge, bad_converge, get_nopsin_eig, get_eigenvalues, final_converge],
                         nupdown_functionals=nupdown_functionals, nupdowns=nupdowns, names=names, vaspobj=vaspobj)

class TSWSBulkChain(SpinCustomChain):
    def __init__(self, vaspobj: Vasp, nupdowns, initial, final, standard=[], override=[], final_step='5_hse' ):
        standard = [load_default_vasp, ws_standard, ws_bulk, load_optimized_U_species, rough_converge, set_222, set_iopt_7, set_single_neb]
        override.append(set_prec_normal)
        gamma = [set_gamma, gamma_optimization]
        pre_converge   = CustomFunctional(Vasp, standard + [awful_converge] + override + gamma, type='neb')
        bad_converge   = CustomFunctional(Vasp, standard + [rough_converge] + override, type='neb')
        bad_converge_gamma   = CustomFunctional(Vasp, standard + [rough_converge] + override + gamma , type='neb')
        get_nopsin_eig = CustomFunctional(Vasp, standard + [get_eigen_nospin] + override, type='neb')
        get_nopsin_eig_gamma = CustomFunctional(Vasp, standard + [get_eigen_nospin] + override + gamma, type='neb')
        get_eigenvalues= CustomFunctional(Vasp, standard + [get_eigen] + override, type='neb')
        get_eigenvalues_gamma = CustomFunctional(Vasp, standard + [get_eigen] + override + gamma, type='neb')
        final_converge = CustomFunctional(Vasp, standard + [full_converge, set_dimer] + override)
        final_converge_gamma = CustomFunctional(Vasp, standard + [full_converge, set_dimer] + override + gamma)
        hse            = CustomFunctional(Vasp, standard + [hse06, set_algo_damp, set_nkred_222] + override + [single_point])
        names          = ['0_pre_converge', '1_rough_converge', '2_nospin_eig', '3_get_eigenvalues', '4_final_converge', final_step]
        nupdown_functionals = [pre_converge, bad_converge_gamma, get_nopsin_eig_gamma, get_eigenvalues_gamma, final_converge_gamma]
        super().__init__([pre_converge, bad_converge, get_nopsin_eig, get_eigenvalues, final_converge, hse], initial_structure=initial.copy(), final_structure=final.copy(),
                          nupdown_functionals=nupdown_functionals, nupdowns=nupdowns, names=names, vaspobj=vaspobj)

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
        dos            = CustomFunctional(Vasp, standard + [single_point, hse06, set_algo_damp, set_dos, tetrahedron, all_output])
        names          = ['0_pre_converge', '1_rough_converge', '2_nospin_eig', '3_get_eigenvalues', '4_final_converge', '5_hse']
        return super().__init__([pre_converge, bad_converge, get_nopsin_eig, get_eigenvalues, final_converge, hse], names=names, vaspobj=vaspobj)

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
        dos            = CustomFunctional(Vasp, standard + [single_point, hse06, set_algo_damp, tetrahedron, all_output])
        names          = ['0_pre_converge', '1_rough_converge', '2_nospin_eig', '3_get_eigenvalues', '4_final_converge', '5_hse']
        return super().__init__([pre_converge, bad_converge, get_nopsin_eig, get_eigenvalues, final_converge, hse], names=names)

class WSSurfaceChain(CustomChain):
    def __init__(self, vaspobj: Vasp, standard = [], override = [], additional_names=[], addition_steps=[]):
        spin = ferro_spin
        standard = [load_default_vasp, ws_standard, ws_surface, load_optimized_U_species, spin, rough_converge, set_221, set_iopt_7, set_kpar_auto] + standard
        pre_converge = CustomFunctional(Vasp, standard + [awful_converge, set_gamma, gamma_optimization, set_algo_fast] + override)
        bad_converge = CustomFunctional(Vasp, standard + [rough_converge, set_algo_fast] + override)
        get_nopsin_eig = CustomFunctional(Vasp, standard + [get_eigen_nospin, set_algo_normal] + override)
        get_eigenvalues = CustomFunctional(Vasp, standard + [get_eigen, set_algo_normal] + override)
        final_converge = CustomFunctional(Vasp, standard + [full_converge, set_algo_fast, all_output ] + override)
#        dos = CustomFunctional(Vasp, standard + [single_point, hse06, set_algo_damp, set_nkred_221, tetrahedron, all_output, set_dos] + override)
        names          = ['0_pre_converge', '1_rough_converge', '2_nospin_eig', '3_get_eigenvalues', '4_final_converge']
        functionals = [pre_converge, bad_converge, get_nopsin_eig, get_eigenvalues, final_converge]
        for i, step in enumerate(addition_steps):
            names += [additional_names[i]]
            functionals += [CustomFunctional(Vasp, standard + step + override)]
        return super().__init__(functionals, names=names, vaspobj=vaspobj)

class WSSurfaceChain_hse(CustomChain):
    def __init__(self, vaspobj: Vasp, standard = [], override = [], additional_names=[], addition_steps=[]):
        spin = ferro_spin
        standard = [load_default_vasp, ws_standard, ws_surface, load_optimized_U_species, spin, rough_converge, set_221, set_iopt_7, set_kpar_auto] + standard
        pre_converge = CustomFunctional(Vasp, standard + [awful_converge, set_gamma, gamma_optimization, set_algo_fast] + override)
        bad_converge = CustomFunctional(Vasp, standard + [rough_converge, set_algo_fast] + override)
        get_nopsin_eig = CustomFunctional(Vasp, standard + [get_eigen_nospin, set_algo_normal] + override)
        get_eigenvalues = CustomFunctional(Vasp, standard + [get_eigen, set_algo_normal] + override)
        final_converge = CustomFunctional(Vasp, standard + [full_converge, set_algo_fast, all_output ] + override)
        hse = CustomFunctional(Vasp, standard + [single_point, hse06, set_algo_damp035, set_nkred_221, all_output] + override)
#        dos = CustomFunctional(Vasp, standard + [single_point, hse06, set_algo_damp, set_nkred_221, tetrahedron, all_output, set_dos] + override)
        names          = ['0_pre_converge', '1_rough_converge', '2_nospin_eig', '3_get_eigenvalues', '4_final_converge', '5_hse']
        functionals    = [pre_converge,      bad_converge,       get_nopsin_eig, get_eigenvalues,     final_converge,     hse]
        for i, step in enumerate(addition_steps):
            names += [additional_names[i]]
            functionals += [CustomFunctional(Vasp, standard + step + override)]
        return super().__init__(functionals, names=names, vaspobj=vaspobj)

class WSSurfaceChain_unit(WSSurfaceChain):
    def __init__(self, vaspobj: Vasp):
        return  super().__init__(vaspobj, standard=[set_kpoints_auto_20])

class WSSurfaceChain_unit_vib(WSSurfaceChain):
    def __init__(self, vaspobj: Vasp):
        additional_steps = [[vibrations_disp, set_algo_fast]]
        return  super().__init__(vaspobj, standard=[set_kpoints_auto_20], addition_steps=additional_steps, additional_names=['5_vibration'])

class WSSurfaceChain_gamma(WSSurfaceChain):
    def __init__(self, vaspobj : Vasp):
        return  super().__init__(vaspobj, standard=[set_gamma, gamma_optimization], override=[unset_nkred])

class WSSurfaceChain_gamma_dimer(WSSurfaceChain):
    def __init__(self, vaspobj : Vasp):
        return  super().__init__(vaspobj, standard=[set_gamma, gamma_optimization, set_dimer], override=[unset_nkred])

class WSBulkPBE(OptimizedParametersChain):
    def __init__(self, vaspobj: Vasp, bandgap:float=None, standard=[], override=[], final_step='5_hse' ):
        standard = [load_default_vasp, cell_relax, ws_bulk, load_optimized_U_species, set_kpar_by_core]
        pbe = CustomFunctional(Vasp, standard)
        pbe_single = CustomFunctional(Vasp, standard + [all_output])
        names = ['1_pbe', '2_pbe_reconverge']
        super().__init__([pbe, pbe_single], bandgap=bandgap, names=names, vaspobj=vaspobj)

class WSBulkSCAN(OptimizedParametersChain):
    def __init__(self, vaspobj: Vasp, bandgap:float=None, standard=[], override=[], final_step='5_hse' ):
        standard = [load_default_vasp, cell_relax, ws_bulk, scan, set_kpar_by_core]
        pbe = CustomFunctional(Vasp, standard)
        pbe_single = CustomFunctional(Vasp, standard + [all_output])
        names = ['1_scan', '2_scan_reconverge']
        super().__init__([pbe, pbe_single], bandgap=bandgap, names=names, vaspobj=vaspobj)

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
    'O'  :  0,
    'Al' : 0,
    'H'  : 0

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
    vasp.encut = 600
    vasp.ediff = 1e-6
    vasp.ispin = 2
    vasp.magmom = True
    vasp.sigma = 0.01
    vasp.ldau = True
    vasp.maxmix = 200
    vasp.isym = 0
    return vasp

def ws_surface(vasp: Vasp, structure):
    ws_bulk(vasp, structure)
    vasp.encut = 800
    return vasp
