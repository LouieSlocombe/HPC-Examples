import ase
from ase.calculators.nwchem import NWChem
from ase.optimize import BFGS


def get_nwchem_calculator(f_charge=0,
                          f_xc_functional="B3LYP",
                          f_multiplicity=1,
                          f_basis_set="6-311++G**",
                          f_disp=None,
                          f_solv=None):
    # Init the nwchem dictionary
    tmp = dict(label='calc/nwchem', charge=f_charge)
    # Add entry for basis set
    tmp["basis"] = f_basis_set

    # Check if using cam-b3lyp
    if f_xc_functional.upper() == "CAM-B3LYP":
        xc_f_cam_b3lyp = "xcamb88 1.00 lyp 0.81 vwn_5 0.19 hfexch 1.00"
        xc_cam = "0.33 cam_alpha 0.19 cam_beta 0.46"
        # Make the dft block
        tmp["dft"] = dict(maxiter=2000,
                          iterations=1000,
                          grid="fine nodisk",
                          direct=" ",
                          noio=" ",
                          xc=xc_f_cam_b3lyp,
                          cam=xc_cam,
                          mult=f_multiplicity)
    else:
        # Make the standard dft block
        tmp["dft"] = dict(maxiter=2000,
                          iterations=1000,
                          grid="fine nodisk",
                          direct=" ",
                          noio=" ",
                          xc=str(f_xc_functional).upper(),
                          mult=f_multiplicity)

    if f_disp.upper() == "XDM":
        # https://nwchemgit.github.io/Density-Functional-Theory-for-Molecules.html#xdm-exchange-hole-dipole-moment-dispersion-model
        val = tmp.get("dft")  # Key the key value
        val["xdm "] = "a1 0.6224 a2 1.7068"  # modify the value
        tmp["dft"] = val  # Put it back

    elif f_disp.upper() == "D3":
        # https://nwchemgit.github.io/Density-Functional-Theory-for-Molecules.html#disp-empirical-long-range-contribution-vdw
        val = tmp.get("dft")  # Key the key value
        val["disp"] = "vdw 3"  # modify the value
        tmp["dft"] = val  # Put it back

    if f_solv.upper() == "WATER":
        # https://nwchemgit.github.io/COSMO-Solvation-Model.html
        tmp["cosmo"] = dict(do_cosmo_smd=True, solvent='water')

    if f_solv.upper() == "PROTEIN":
        # https://nwchemgit.github.io/COSMO-Solvation-Model.html
        tmp["cosmo"] = dict(do_cosmo_smd=True, dielec=8.0)

    return NWChem(**tmp)


calc = get_nwchem_calculator(f_charge=1,
                             f_xc_functional="B3LYP",
                             f_multiplicity=1,
                             f_basis_set="6-311++G**",
                             f_disp="D3",
                             f_solv="PROTEIN")
# load the structure
# atoms = ase.io.mol.read_mol("water.mol")
atoms = ase.io.read("NaH2O4.traj")
# set the calculator
atoms.calc = calc
# optimise the structure
Optimise = ase.optimize.BFGS(atoms, trajectory="opti.traj", logfile="opti.log")
Optimise.run(fmax=0.01)
