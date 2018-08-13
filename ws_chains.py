from default_chains import *
from pylada.vasp.specie import U

def load_optimized_U_species(vasp : Vasp, structure):
    # See vasp/functional.py:  elementName, fileName, max or min oxidation state
    pseudoDir = '$PSEUDO_DIR'
    vasp.add_specie = "Sc", pseudoDir + "/Sc_sv", U("dudarev", "d", 3)    # FERE
    vasp.add_specie = "Y", pseudoDir + "/Y_sv", U("dudarev", "d", 3)    # FERE
    vasp.add_specie = "Ti", pseudoDir + "/Ti_pv", U("dudarev", "d", 4.35) # Wolverton 3+
    vasp.add_specie = "V",  pseudoDir + "/V_pv",  U("dudarev", "d", 4.9)  # Wolverton 2+
    vasp.add_specie = "Cr", pseudoDir + "/Cr_pv", U("dudarev", "d", 3.04) # Wolverton 3+
    vasp.add_specie = "Zr", pseudoDir + "/Zr_sv", U("dudarev", "d", 3) # FERE
    vasp.add_specie = "Mn", pseudoDir + "/Mn_pv", U("dudarev", "d", 2.98) # Wolverton 2+
    vasp.add_specie = "Mn3p", pseudoDir + "/Mn_pv", U("dudarev", "d", 4.54) # Wolverton 3+
    vasp.add_specie = "Fe", pseudoDir + "/Fe_pv", U("dudarev", "d", 4.04) # Wolverton 2+,
    vasp.add_specie = "Co", pseudoDir + "/Co_pv", U("dudarev", "d", 3.75) # Wolverton 2+
    vasp.add_specie = "Ni",  pseudoDir + "/Ni_pv", U("dudarev", "d", 4.4) # Wolverton 2+
    vasp.add_specie = "Cu", pseudoDir + "/Cu_pv", U("dudarev", "d" , 5.0) # FERE
    vasp.add_specie = "Nb", pseudoDir + "/Nb_pv", U("dudarev", "d" , 3.0) # FERE
    vasp.add_specie = "Ce", pseudoDir + "/Ce", U("dudarev", "f" , 3.0) # Wolverton Ceria Paper
    vasp.add_specie = "Zn", pseudoDir + "/Zn",
    vasp.add_specie = "Ga", pseudoDir + "/Ga",
    vasp.add_specie = "O",  pseudoDir + "/O"
    vasp.add_specie = "Al", pseudoDir + "/Al"
    vasp.add_specie = "H", pseudoDir + "/H"
    
    vasp.add_specie = "Sr", pseudoDir + "/Sr_sv"
    vasp.add_specie = "Si", pseudoDir + "/Si"
    vasp.add_specie = "Cd", pseudoDir + "/Cd"
    vasp.add_specie = "Ba", pseudoDir + "/Ba_sv"
    vasp.add_specie = "Be", pseudoDir + "/Be"
    vasp.add_specie = "Ba", pseudoDir + "/Ba_sv"
    vasp.add_specie = "Bi", pseudoDir + "/Bi_d"
    vasp.add_specie = "Ge", pseudoDir + "/Ge_d"
    vasp.add_specie = "Se", pseudoDir + "/Se"
    vasp.add_specie = "Sb", pseudoDir + "/Sb"
    vasp.add_specie = "Na", pseudoDir + "/Na_pv"
    vasp.add_specie = "Ca", pseudoDir + "/Ca_pv"
    vasp.add_specie = "In", pseudoDir + "/In_d"
    vasp.add_specie = "K", pseudoDir + "/K_sv"
    vasp.add_specie = "Mg", pseudoDir + "/Mg"
    vasp.add_specie = "Sn", pseudoDir + "/Sn_d"
    vasp.add_specie = "Li", pseudoDir + "/Li_sv"

    vasp.add_specie = "La", pseudoDir + "/La" # TODO Determine U
    return(vasp)

def load_low_FERE_species(vasp: Vasp, structure):
    pseudoDir = '$PSEUDO_DIR'
    vasp.add_specie = "Fe", pseudoDir + "/Fe", U("dudarev", "d", 3)  # FERE
    vasp.add_specie = "O", pseudoDir + "/O_s"
    vasp.add_specie = "Al", pseudoDir + "/Al"
    return vasp


def load_FERE_species(vasp: Vasp, structure):
    # See vasp/functional.py:  elementName, fileName, max or min oxidation state
    pseudoDir = '$PSEUDO_DIR'
    vasp.add_specie = "Sc", pseudoDir + "/Sc_sv", U("dudarev", "d", 3)  # FERE
    vasp.add_specie = "Y", pseudoDir + "/Y_sv", U("dudarev", "d", 3)  # FERE
    # vasp.add_specie = "Ti", pseudoDir + "/Ti_pv", U("dudarev", "d", 4.35)  # Wolverton 3+
    # vasp.add_specie = "V", pseudoDir + "/V_pv", U("dudarev", "d", 4.9)  # Wolverton 2+
    # vasp.add_specie = "Cr", pseudoDir + "/Cr_pv", U("dudarev", "d", 3.04)  # Wolverton 3+
    vasp.add_specie = "Zr", pseudoDir + "/Zr_sv", U("dudarev", "d", 3)  # FERE
    # vasp.add_specie = "Mn", pseudoDir + "/Mn_pv", U("dudarev", "d", 2.98)  # Wolverton 2+
    # vasp.add_specie = "Mn3p", pseudoDir + "/Mn_pv", U("dudarev", "d", 4.54)  # Wolverton 3+
    vasp.add_specie = "Fe", pseudoDir + "/Fe_pv", U("dudarev", "d", 3)  # FERE
    vasp.add_specie = "Co", pseudoDir + "/Co_pv", U("dudarev", "d", 3)  # Wolverton 2+
    vasp.add_specie = "Ni", pseudoDir + "/Ni_pv", U("dudarev", "d", 3)  # FERE
    vasp.add_specie = "Cu", pseudoDir + "/Cu_pv", U("dudarev", "d", 5.0)  # FERE
    vasp.add_specie = "Nb", pseudoDir + "/Nb_pv", U("dudarev", "d", 3.0)  # FERE
    vasp.add_specie = "Ce", pseudoDir + "/Ce", U("dudarev", "f", 3.0)  # Wolverton Ceria Paper
    vasp.add_specie = "Zn", pseudoDir + "/Zn",
    vasp.add_specie = "Ga", pseudoDir + "/Ga",
    vasp.add_specie = "O", pseudoDir + "/O"
    vasp.add_specie = "Al", pseudoDir + "/Al"
    vasp.add_specie = "H", pseudoDir + "/H"

    vasp.add_specie = "Sr", pseudoDir + "/Sr_sv"
    vasp.add_specie = "Si", pseudoDir + "/Si"
    vasp.add_specie = "Cd", pseudoDir + "/Cd"
    vasp.add_specie = "Ba", pseudoDir + "/Ba_sv"
    vasp.add_specie = "Be", pseudoDir + "/Be"
    vasp.add_specie = "Ba", pseudoDir + "/Ba_sv"
    vasp.add_specie = "Bi", pseudoDir + "/Bi_d"
    vasp.add_specie = "Ge", pseudoDir + "/Ge_d"
    vasp.add_specie = "Se", pseudoDir + "/Se"
    vasp.add_specie = "Sb", pseudoDir + "/Sb"
    vasp.add_specie = "Na", pseudoDir + "/Na_pv"
    vasp.add_specie = "Ca", pseudoDir + "/Ca_pv"
    vasp.add_specie = "In", pseudoDir + "/In_d"
    vasp.add_specie = "K", pseudoDir + "/K_sv"
    vasp.add_specie = "Mg", pseudoDir + "/Mg"
    vasp.add_specie = "Sn", pseudoDir + "/Sn_d"
    vasp.add_specie = "Li", pseudoDir + "/Li_sv"

    vasp.add_specie = "La", pseudoDir + "/La"  # TODO Determine U
    return (vasp)

class WSBulkChain(SpinCustomChain):
    def __init__(self, vaspobj: Vasp(), nupdowns=[], standard=[], override=[], final_step='5_single_point' ):
        standard = [load_default_vasp, ws_standard, herc_bulk, load_optimized_U_species, rough_converge, set_222, set_iopt_7]
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
        hse            = CustomFunctional(Vasp, standard + [full_converge] + override + [single_point])
        names          = ['0_pre_converge', '1_rough_converge', '2_nospin_eig', '3_get_eigenvalues', '4_final_converge', final_step]
        nupdown_functionals = [pre_converge, bad_converge_gamma, get_nopsin_eig_gamma, get_eigenvalues_gamma, final_converge_gamma]
        super().__init__([pre_converge, bad_converge, get_nopsin_eig, get_eigenvalues, final_converge, hse],
                         nupdown_functionals=nupdown_functionals, nupdowns=nupdowns, names=names, vaspobj=vaspobj)

class WSBulkChain_FERE(WSBulkChain):
    def __init__(self, vaspobj: Vasp(), nupdowns=[], standard=[], override=[], final_step='5_single_point'):
        super().__init__(vaspobj, nupdowns=nupdowns, standard=standard, override=[load_low_FERE_species, anti_spin_interstitial], final_step=final_step)

class WSBulkChain_auto(SpinCustomChain):
    def __init__(self, vaspobj: Vasp(), nupdowns, standard=[], override=[], **kwargs):
        standard = [load_default_vasp, ws_standard, ws_bulk, load_optimized_U_species, rough_converge, set_iopt_7, set_kpar_auto, set_spin]
        gamma = [set_gamma, gamma_optimization]
        pre_converge   = CustomFunctional(Vasp, standard + [awful_converge, set_algo_fast] + override + gamma)
        bad_converge   = CustomFunctional(Vasp, standard + [rough_converge, set_algo_fast] + override)
        bad_converge_gamma   = CustomFunctional(Vasp, standard + [rough_converge, set_algo_fast] + gamma + override)
        get_nopsin_eig = CustomFunctional(Vasp, standard + [get_eigen_nospin] + override)
        get_nopsin_eig_gamma = CustomFunctional(Vasp, standard + [get_eigen_nospin] + gamma + override)
        get_eigenvalues= CustomFunctional(Vasp, standard + [get_eigen] + override)
        get_eigenvalues_gamma = CustomFunctional(Vasp, standard + [get_eigen] + gamma + override)
        final_converge = CustomFunctional(Vasp, standard + [full_converge,all_output] + override)
        final_converge_gamma = CustomFunctional(Vasp, standard + [full_converge] + gamma + override)
        sp = CustomFunctional(Vasp, standard + [full_converge, all_output] + override)

        names          = ['0_pre_converge', '1_rough_converge', '2_nospin_eig', '3_get_eigenvalues', '4_final_converge', '5_single_point']
        nupdown_functionals = [pre_converge, bad_converge_gamma, get_nopsin_eig_gamma, get_eigenvalues_gamma, final_converge_gamma]
        super().__init__([pre_converge, bad_converge, get_nopsin_eig, get_eigenvalues, final_converge, sp],
                         nupdown_functionals=nupdown_functionals, nupdowns=nupdowns, names=names, vaspobj=vaspobj, **kwargs)

class WSBulkChainNoNUPDOWN(CustomChain):
    def __init__(self, vaspobj=Vasp(), standard=[], override=[], **kwargs):
        standard = [load_default_vasp, ws_standard, ws_surface, load_optimized_U_species, rough_converge, set_iopt_7, set_kpar_auto, set_kpoints_auto_20]
        gamma = [set_gamma, gamma_optimization]
        pre_converge   = CustomFunctional(Vasp, standard + [awful_converge, set_algo_fast] + gamma + override)
        bad_converge   = CustomFunctional(Vasp, standard + [rough_converge, set_algo_fast] + override)
        get_nopsin_eig = CustomFunctional(Vasp, standard + [get_eigen_nospin] + override)
        get_eigenvalues= CustomFunctional(Vasp, standard + [get_eigen] + override)
        final_converge = CustomFunctional(Vasp, standard + [full_converge,all_output] + override)
        sp = CustomFunctional(Vasp, standard + [full_converge, all_output] + override)

        names          = ['0_pre_converge', '1_rough_converge', '2_nospin_eig', '3_get_eigenvalues', '4_final_converge', '5_single_point']
        super().__init__([pre_converge, bad_converge, get_nopsin_eig, get_eigenvalues, final_converge, sp],
                         names=names, vaspobj=vaspobj, **kwargs)

class WSSurfaceIntermediateChain(SpinCustomChain):
    def __init__(self, vaspobj: Vasp(), nupdowns, standard=[], override=[], **kwargs):
        standard = [load_default_vasp, ws_standard, ws_surface, load_optimized_U_species, rough_converge, set_iopt_7, set_kpar_auto, set_kpoints_auto_20]
        gamma = [set_gamma, gamma_optimization]
        pre_converge   = CustomFunctional(Vasp, standard + [awful_converge, set_algo_fast] + gamma + override)
        bad_converge   = CustomFunctional(Vasp, standard + [rough_converge, set_algo_fast] + override)
        bad_converge_gamma   = CustomFunctional(Vasp, standard + [rough_converge, set_algo_fast] + gamma + override)
        get_nopsin_eig = CustomFunctional(Vasp, standard + [get_eigen_nospin] + override)
        get_nopsin_eig_gamma = CustomFunctional(Vasp, standard + [get_eigen_nospin] + gamma + override)
        get_eigenvalues= CustomFunctional(Vasp, standard + [get_eigen] + override)
        get_eigenvalues_gamma = CustomFunctional(Vasp, standard + [get_eigen] + gamma + override)
        final_converge = CustomFunctional(Vasp, standard + [full_converge,all_output,set_ediffg_03] + override)
        final_converge_gamma = CustomFunctional(Vasp, standard + [full_converge, set_low_ediffg] + gamma + override)
        sp = CustomFunctional(Vasp, standard + [full_converge, all_output] + override)

        names          = ['0_pre_converge', '1_rough_converge', '2_nospin_eig', '3_get_eigenvalues', '4_final_converge', '5_single_point']
        nupdown_functionals = [pre_converge, bad_converge_gamma, get_nopsin_eig_gamma, get_eigenvalues_gamma, final_converge_gamma]
        super().__init__([pre_converge, bad_converge, get_nopsin_eig, get_eigenvalues, final_converge, sp],
                         nupdown_functionals=nupdown_functionals, nupdowns=nupdowns, names=names, vaspobj=vaspobj, **kwargs)

class WSSurfaceIntermediateChain_low(SpinCustomChain):
    def __init__(self, vaspobj: Vasp(), nupdowns, standard=[], override=[], **kwargs):
        standard = [load_default_vasp, ws_standard, ws_surface, load_optimized_U_species, rough_converge, set_iopt_7, set_kpar_auto, set_kpoints_auto_20]
        gamma = [set_gamma, gamma_optimization]
        pre_converge   = CustomFunctional(Vasp, standard + [awful_converge, set_algo_fast] + gamma + override)
        bad_converge_gamma   = CustomFunctional(Vasp, standard + [rough_converge, set_algo_fast] + gamma + override)
        get_nopsin_eig_gamma = CustomFunctional(Vasp, standard + [get_eigen_nospin] + gamma + override)
        get_eigenvalues_gamma = CustomFunctional(Vasp, standard + [get_eigen, set_prec_normal] + gamma + override)
        final_converge_gamma = CustomFunctional(Vasp, standard + [full_converge, set_low_ediffg, set_prec_normal] + gamma + override)

        names          = ['lowest_nupdown']
        nupdown_functionals = [pre_converge, bad_converge_gamma, get_nopsin_eig_gamma, get_eigenvalues_gamma, final_converge_gamma]
        super().__init__([final_converge_gamma],
                         nupdown_functionals=nupdown_functionals, nupdowns=nupdowns, names=names, vaspobj=vaspobj, **kwargs)

class WSBulkChain_small(SpinCustomChain):
    def __init__(self, vaspobj: Vasp, nupdowns, standard=[], override=[], final_step='5_hse'):
        standard = [load_default_vasp, ws_standard, herc_bulk, load_optimized_U_species, rough_converge, set_444, set_iopt_7]
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
        standard = [load_default_vasp, ws_standard, herc_bulk, load_optimized_U_species, rough_converge, set_222, set_iopt_7, set_single_neb]
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
        standard = [load_default_vasp, ws_standard, herc_bulk, load_optimized_U_species, spin, rough_converge, set_222, set_iopt_7]
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
        standard = [load_default_vasp, ws_standard, herc_bulk, load_optimized_U_species, spin, rough_converge, set_222, set_iopt_7]
        pre_converge   = CustomFunctional(Vasp, standard + [awful_converge, set_gamma, gamma_optimization])
        bad_converge   = CustomFunctional(Vasp, standard + [rough_converge])
        get_nopsin_eig = CustomFunctional(Vasp, standard + [get_eigen_nospin, set_algo_normal])
        get_eigenvalues= CustomFunctional(Vasp, standard + [get_eigen])
        final_converge = CustomFunctional(Vasp, standard + [full_converge])
        hse            = CustomFunctional(Vasp, standard + [hse06, set_nkred_222])
        dos            = CustomFunctional(Vasp, standard + [single_point, hse06, set_algo_damp, tetrahedron, all_output])
        names          = ['0_pre_converge', '1_rough_converge', '2_nospin_eig', '3_get_eigenvalues', '4_final_converge', '5_hse']
        return super().__init__([pre_converge, bad_converge, get_nopsin_eig, get_eigenvalues, final_converge, hse], names=names)

class WSSurfaceChain(SpinCustomChain):
    def __init__(self, vaspobj: Vasp, nupdowns, standard=[], override=[], additional_names=[], addition_steps=[]):
        spin = ferro_spin
        gamma = [set_gamma, gamma_optimization]
        standard = [load_default_vasp, ws_standard, herc_surface, load_optimized_U_species, spin, rough_converge, set_221, set_iopt_7, set_kpar_auto] + standard
        pre_converge = CustomFunctional(Vasp, standard + [awful_converge, set_gamma, gamma_optimization, set_algo_fast] + override)
        bad_converge = CustomFunctional(Vasp, standard + [rough_converge, set_algo_fast] + override)
        get_nopsin_eig = CustomFunctional(Vasp, standard + [get_eigen_nospin, set_algo_normal] + override)
        get_eigenvalues = CustomFunctional(Vasp, standard + [get_eigen, set_algo_normal] + override)
        final_converge = CustomFunctional(Vasp, standard + [full_converge, set_algo_fast, all_output ] + override)
        dos = CustomFunctional(Vasp, standard + [single_point, full_converge, set_algo_normal, set_nkred_221, tetrahedron, all_output, set_dos] + override)
        names          = ['0_pre_converge', '1_rough_converge', '2_nospin_eig', '3_get_eigenvalues', '4_final_converge', '5_dos']
        functionals = [pre_converge, bad_converge, get_nopsin_eig, get_eigenvalues, final_converge, dos]
        for i, step in enumerate(addition_steps):
            names += [additional_names[i]]
            functionals += [CustomFunctional(Vasp, standard + step + override)]

        bad_converge_gamma = CustomFunctional(Vasp, standard + [rough_converge, set_algo_fast] + override + gamma)
        get_nopsin_eig_gamma = CustomFunctional(Vasp, standard + [get_eigen_nospin, set_algo_normal] + override + gamma)
        get_eigenvalues_gamma = CustomFunctional(Vasp, standard + [get_eigen, set_algo_normal] + override + gamma)
        final_converge_gamma = CustomFunctional(Vasp, standard + [full_converge, set_algo_fast, all_output] + override + gamma)
        nupdown_functionals = [pre_converge, bad_converge_gamma, get_nopsin_eig_gamma, get_eigenvalues_gamma, final_converge_gamma]
        return super().__init__(functionals, nupdown_functionals=nupdown_functionals, nupdowns=nupdowns, names=names, vaspobj=vaspobj)

class WSSurfaceChain_hse(CustomChain):
    def __init__(self, vaspobj: Vasp, standard = [], override = [], additional_names=[], addition_steps=[]):
        spin = ferro_spin
        standard = [load_default_vasp, ws_standard, herc_surface, load_optimized_U_species, spin, rough_converge, set_221, set_iopt_7, set_kpar_auto] + standard
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
        standard = [load_default_vasp, cell_relax, herc_bulk, load_optimized_U_species, set_kpar_by_core, set_spin]
        pbe = CustomFunctional(Vasp, standard)
        pbe_single = CustomFunctional(Vasp, standard + [all_output])
        names = ['1_pbe', '2_pbe_reconverge']
        super().__init__([pbe, pbe_single], bandgap=bandgap, names=names, vaspobj=vaspobj)

class WSBulkToSurfacePBE(CustomChain):
    def __init__(self, vaspobj: Vasp, bulk_structure, standard=[], override=[], incar_settings='../INCAR.defaults', kpt_modifier=1):
        from Helpers import pyl_to_pmg
        standard = [load_default_vasp, ws_bulk, load_optimized_U_species, set_kpar_2, set_iopt_7, idipol_3, set_algo_normal, set_nelm_200, set_iopt_7]
        with open(incar_settings) as f:
            lines = [line.strip().split('=') for line in f.readlines()]
            incar = {f[0].strip(): float(f[1]) for f in lines}
            kpts = math.ceil((incar['KPOINTS'] - 0.25) * max(pyl_to_pmg(bulk_structure).lattice.abc)/ kpt_modifier)

        pre_converge = CustomFunctional(Vasp, standard + [awful_converge, set_gamma, gamma_optimization] + override)
        bad_converge = CustomFunctional(Vasp, standard + [rough_converge] + override)
        get_nopsin_eig = CustomFunctional(Vasp, standard + [get_eigen_nospin, no_relax] + override)
        get_eigenvalues = CustomFunctional(Vasp, standard + [get_eigen] + override)
        final_converge = CustomFunctional(Vasp, standard + [full_converge, all_output] + override)
        ldipol = CustomFunctional(Vasp, standard + [full_converge, all_output, surface_final] + override)
        names = ['0_pre_converge', '1_rough_converge', '2_nospin_eig', '3_get_eigenvalues', '4_final_converge', '5_ldipol']
        functionals = [pre_converge, bad_converge, get_nopsin_eig, get_eigenvalues, final_converge, ldipol]
        super().__init__(functionals, names=names, vaspobj=vaspobj, encut=incar['ENCUT'], kpoints=kpts)

class WSBulkToFrozenSurfacePBE(CustomChain):
    def __init__(self, vaspobj: Vasp, bulk_structure, standard=[], override=[], incar_settings='../INCAR.defaults', kpt_modifier=1):
        from Helpers import pyl_to_pmg
        standard = [load_default_vasp, ws_bulk, load_optimized_U_species, set_kpar_2, idipol_3, no_relax]
        with open(incar_settings) as f:
            lines = [line.strip().split('=') for line in f.readlines()]
            incar = {f[0].strip(): float(f[1]) for f in lines}
            kpts = math.ceil((incar['KPOINTS'] - 0.25) * max(pyl_to_pmg(bulk_structure).lattice.abc) / kpt_modifier)

        get_nopsin_eig = CustomFunctional(Vasp, standard + [get_eigen_nospin, set_algo_fast, set_nelm_200] + override)
        get_eigenvalues = CustomFunctional(Vasp, standard + [get_eigen, set_nelm_200] + override)
        final_converge = CustomFunctional(Vasp, standard + [full_converge, set_algo_normal, all_output, set_nelm_200] + override)
        ldipol = CustomFunctional(Vasp, standard + [full_converge, set_algo_fast, all_output, set_nelm_200, surface_final] + override)
        names = ['2_nospin_eig', '3_get_eigenvalues', '4_final_converge', '5_ldipol']
        functionals = [get_nopsin_eig, get_eigenvalues, final_converge, ldipol]
        super().__init__(functionals, names=names, vaspobj=vaspobj, encut=incar['ENCUT'], kpoints=kpts)

def make_surfaces_to_pylada(root, bulk_structure, incar_settings=None, label='', depth=8, frozen_depth=2, override=[], kpt_modifier=1):
    '''

    :param root: root directory of run
    :param bulk_structure: pylada bulk strucutre
    :param incar_settings: location of INCAR.defaults
    :return:
    '''
    from pylada.crystal import write
    from Generate_Surface import Generate_Surface
    from Helpers import pyl_to_pmg, pmg_to_pyl
    from Generate_Surface import get_bottom, get_SD_along_vector
    # small_surfaces = Generate_Surface(pyl_to_pmg(bulk_structure), 1, 1, 1, 3, vacuum=8, orth=True)
    for i, surface in enumerate(Generate_Surface(pyl_to_pmg(bulk_structure), 1, 1, 1, depth, vacuum=8, orth=False)):
        # Frozen Surface
        # surface_small = small_surfaces[i]
        surf_folder = root / label / str(i).zfill(2)
        surf_folder.functional = WSBulkToFrozenSurfacePBE(Vasp(), bulk_structure=bulk_structure,
                                                          incar_settings=incar_settings, override=override,
                                                          kpt_modifier=kpt_modifier)
        surf_folder.params['structure'] = pmg_to_pyl(surface).copy()
        os.makedirs(surf_folder.name[1:], exist_ok=True)
        surface.to('poscar', os.path.join(surf_folder.name[1:], 'surface.vasp'))
        with open(os.path.join(surf_folder.name[1:], 'DATABASE'), 'w') as f:
            f.write('''
surface
surface_cut {}
surface_termination  {}
convergence_study
convergence_type surface
misc_labels {}
'''.format(i, 'None', label))

        # Frozen Surfaces

        for frozen_region in ['top', 'bottom']:

            froz_folder = surf_folder / frozen_region
            froz_folder.functional = WSBulkToSurfacePBE(Vasp(), bulk_structure=bulk_structure,
                                                        incar_settings=incar_settings, override=override,
                                                        kpt_modifier=kpt_modifier)

            surface_frozen = surface.copy()
            # print(surface_frozen)
            sd = get_SD_along_vector(surface_frozen, 2, get_bottom(surface_frozen, length=frozen_depth, region=frozen_region))
            # Poscar(surface_frozen, selective_dynamics=sd).write_file(os.path.join(froz_folder.name[1:], 'surface.vasp'))
            surface_frozen_pyl = pmg_to_pyl(surface_frozen)
            # sd_small = get_SD_along_vector(surface_frozen, frozen_depth, get_bottom(surface_frozen, region=frozen_region))
            for (atom, sd) in zip(surface_frozen_pyl, sd):
                if sd[0]:
                    atom.freeze = 'xyz'
            os.makedirs(froz_folder.name[1:], exist_ok=True)
            write.poscar(surface_frozen_pyl, os.path.join(froz_folder.name[1:], 'surface.vasp'), vasp5=True)
            froz_folder.params['structure'] = surface_frozen_pyl.copy()
            os.makedirs(froz_folder.name[1:], exist_ok=True)
            with open(os.path.join(froz_folder.name[1:], 'DATABASE'), 'w') as f:
                f.write('''
surface
surface_cut {}
surface_termination  {}
misc_labels {}
'''.format(i, frozen_region, label))



class WSBulkSCAN(OptimizedParametersChain):
    def __init__(self, vaspobj: Vasp, bandgap:float=None, standard=[], override=[], final_step='5_hse' ):
        standard = [load_default_vasp, cell_relax, herc_bulk, scan, set_kpar_by_core]
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
    'Mn3p' : 5,
    'Fe' : 4,
    'Co' : 3,
    'Ni' : 2,
    'Cu' : 1,
    'Zn' : 0,
    'O'  :  0,
    'Al' : 0,
    'H'  : 0
}

def set_spin(vasp: Vasp, structure):
    magmom = []
    for atom in structure:
        try:
            magmom.append(round(atom.magmom,1))
        except:
            if atom.type in spins:
                try:
                    magmom.append(spins[atom.type]*atom.spin_dir)
                except AttributeError:
                    magmom.append(spins[atom.type])
            else:
                magmom.append(0)
    magmom_shortened = ''
    prev_spin = magmom[0]
    num_spins = 0
    for spin in magmom:
        if spin == prev_spin:
            num_spins += 1
        else:
            magmom_shortened += '{}*{} '.format(num_spins, prev_spin)
            prev_spin = spin
            num_spins = 1
    magmom_shortened += '{}*{} '.format(num_spins, prev_spin)

    vasp.ispin=2
    vasp.magmom = magmom_shortened
    return vasp

def ws_standard(vasp: Vasp, structure):
    vasp.isym = 0
    vasp.ismear = -5
    return vasp

def anti_spin(vasp, structure):
    vasp.ispin = 2
    vasp.magmom = '32*0 8*4 8*-4 64*0'
    return vasp

def anti_spin_interstitial(vasp, structure):
    vasp.ispin = 2
    vasp.magmom = '32*0 8*4 8*-4 65*0'
    return vasp

def ferro_spin(vasp, structure):
    vasp.ispin = 2
    vasp.magmom = True
    for a in structure:
        if a.type in spins:
            a.magmom = spins[a.type]
    return vasp

def herc_bulk(vasp: Vasp, structure):
    vasp.prec = 'Accurate'
    vasp.encut =  600
    vasp.ediff = 1e-6
    vasp.ispin = 2
    vasp.magmom = True
    vasp.sigma = 0.01
    vasp.ldau = True
    vasp.maxmix = 200
    vasp.isym = 0
    return vasp

def herc_surface(vasp: Vasp, structure):
    herc_bulk(vasp, structure)
    vasp.encut = 800
    return vasp

def ws_bulk(vasp: Vasp, structure):
    vasp.prec = 'Accurate'
    vasp.ediff = 1e-6
    vasp.ispin = 2
    vasp.magmom = True
    vasp.sigma = 0.01
    vasp.ldau = True
    vasp.maxmix = 200
    vasp.isym = 0
    return vasp

def ws_surface(vasp: Vasp, structure):
    vasp.prec = 'Accurate'
    vasp.ediff = 1e-6
    vasp.ispin = 2
    vasp.magmom = True
    vasp.sigma = 0.01
    vasp.ldau = True
    vasp.maxmix = 200
    vasp.isym = 0
    vasp.encut = 600
    return vasp