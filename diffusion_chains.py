from default_chains import *
from ws_chains import load_optimized_species_no_U

class HDiffusionSCANChain(SpinCustomChain):
    def __init__(self, vaspobj: Vasp(), nupdowns, standard=[], override=[], **kwargs):
        standard = [load_default_vasp, diffusion_standard, load_optimized_species_no_U, rough_converge, set_iopt_7, set_kpar_auto, scan]
        gamma = [set_gamma, gamma_optimization, set_ncore_auto]
        pre_converge   = CustomFunctional(Vasp, standard + [awful_converge, set_algo_fast] + gamma + override)
        bad_converge   = CustomFunctional(Vasp, standard + [rough_converge, set_algo_all] + override)
        bad_converge_gamma   = CustomFunctional(Vasp, standard + [rough_converge, set_algo_all] + gamma + override)
        get_nopsin_eig = CustomFunctional(Vasp, standard + [get_eigen_nospin, set_algo_all] + override)
        get_nopsin_eig_gamma = CustomFunctional(Vasp, standard + [get_eigen_nospin, set_algo_all] + gamma + override)
        get_eigenvalues= CustomFunctional(Vasp, standard + [get_eigen, set_algo_all] + override)
        get_eigenvalues_gamma = CustomFunctional(Vasp, standard + [get_eigen, set_algo_all] + gamma + override)
        final_converge = CustomFunctional(Vasp, standard + [full_converge, all_output] + override)
        final_converge_gamma = CustomFunctional(Vasp, standard + [full_converge] + gamma + override)
        sp = CustomFunctional(Vasp, standard + [full_converge, all_output] + override)

        names          = ['0_pre_converge', '1_rough_converge', '2_nospin_eig', '3_get_eigenvalues', '4_final_converge', '5_single_point']
        nupdown_functionals = [pre_converge, bad_converge_gamma, get_nopsin_eig_gamma, get_eigenvalues_gamma, final_converge_gamma]
        super().__init__([pre_converge, bad_converge, get_nopsin_eig, get_eigenvalues, final_converge, sp],
                         nupdown_functionals=nupdown_functionals, nupdowns=nupdowns, names=names, vaspobj=vaspobj, **kwargs)

class NaDiffusionUnitSCANChain(SpinCustomChain):
    def __init__(self, vaspobj: Vasp(), standard=[], override=[], **kwargs):
        standard = [load_default_vasp, diffusion_standard, load_optimized_species_no_U, rough_converge, set_iopt_7, set_kpar_auto, scan, set_npar_2, set_kpar_2, cell_relax]
        gamma = [set_gamma, gamma_optimization, set_ncore_auto]
        get_eigenvalues= CustomFunctional(Vasp, standard + [get_eigen, set_algo_all, dont_continue] + override)
        get_eigenvalues_gamma = CustomFunctional(Vasp, standard + [get_eigen, set_algo_all, dont_continue] + gamma + override)
        final_converge = CustomFunctional(Vasp, standard + [full_converge, all_output] + override)
        final_converge_gamma = CustomFunctional(Vasp, standard + [full_converge, all_output] + gamma + override)

        names          = ['0_converge', '1_single_point']
        nupdown_functionals = [get_eigenvalues_gamma, final_converge_gamma]
        super().__init__([get_eigenvalues, final_converge],
                         nupdown_functionals=nupdown_functionals, names=names, vaspobj=vaspobj, **kwargs)

class HDiffusionUnitSCANChain(SpinCustomChain):
    def __init__(self, vaspobj: Vasp(), standard=[], override=[], **kwargs):
        standard = [load_default_vasp, diffusion_standard, load_optimized_species_no_U, rough_converge, set_iopt_7, set_kpar_auto, scan, set_npar_2, set_kpar_2, cell_relax]
        gamma = [set_gamma, gamma_optimization, set_ncore_auto]
        get_eigenvalues= CustomFunctional(Vasp, standard + [get_eigen, set_algo_all, dont_continue] + override)
        get_eigenvalues_gamma = CustomFunctional(Vasp, standard + [get_eigen, set_algo_all, dont_continue] + gamma + override)
        final_converge = CustomFunctional(Vasp, standard + [full_converge, all_output] + override)
        final_converge_gamma = CustomFunctional(Vasp, standard + [full_converge, all_output] + gamma + override)

        names          = ['0_converge', '1_single_point']
        nupdown_functionals = [get_eigenvalues_gamma, final_converge_gamma]
        super().__init__([get_eigenvalues, final_converge],
                         nupdown_functionals=nupdown_functionals, names=names, vaspobj=vaspobj, **kwargs)

def diffusion_standard(vasp, structure=None):
    vasp.ediff = 1e-5
    vasp.encut = 520
    vasp.ismear = 0
    vasp.sigma = 0.1
    vasp.add_keyword('METAGGA', 'Scan')
    vasp.prec = 'Accurate'
    vasp.ispin = 2
    vasp.magmom = True
    vasp.maxmix = 200
    vasp.isym = 0
    vasp = set_kpoints_1704(vasp, structure)
    return vasp