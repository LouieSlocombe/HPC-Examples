import copy
import os
import time

import ase.calculators.castep
from ase import Atoms
from ase.calculators.emt import EMT
from ase.io import read, write
from ase.lattice.hexagonal import Graphene
from ase.visualize import view


def castep_calc(atoms, xc_f="PBE", e_cut=800, kpts=[4, 4, 1], directory='data', ps_type='ext', ps_ext='.usp'):
    calc = ase.calculators.castep.Castep(keyword_tolerance=1)

    # include interface settings in .param file
    calc._export_settings = True
    calc._pedantic = True

    # reuse the same directory
    calc._seed = directory
    calc._label = directory
    calc._directory = directory
    calc._rename_existing_dir = False

    # calc._link_pspots = True
    # calc._copy_pspots = False
    # calc._build_missing_pspots = False

    # Param settings
    calc.param.reuse = True
    calc.param.xc_functional = xc_f
    calc.param.cut_off_energy = e_cut
    calc.param.finite_basis_corr = 0
    calc.param.num_dump_cycles = 0  # Prevent CASTEP from writing *wvfn* files
    calc.param.write_checkpoint = 'minimal'
    calc.param.max_scf_cycles = 100
    calc.param.opt_strategy = "speed"
    calc.param.mixing_scheme = "Broyden"  # Broyden Pulay
    calc.param.spin_polarised = False

    calc.param.sedc_apply = True
    calc.param.sedc_scheme = 'MBD*'
    calc.param.perc_extra_bands = 50
    calc.param.elec_energy_tol = 1e-6

    # Cell settings
    calc.cell.kpoint_mp_grid = kpts
    calc.cell.fix_com = False
    calc.cell.symmetry_tol = 0.001
    calc.cell.fix_all_cell = True
    calc.cell.symmetry_generate = True
    calc.cell.symmetry_ops = "1.0 0.0 0.0"

    ele_list = atoms.get_chemical_symbols()
    tmp = None
    ps_suff = None
    if xc_f == 'PBE' or xc_f == 'PBE0':
        tmp = 'PBE'
    elif xc_f == 'B3LYP' or xc_f == 'BLYP':
        tmp = 'BLYP'
    elif xc_f == 'LDA':
        tmp = 'LDA'
    else:
        exit('problem with xc and not supporting ps type...')
    ps_type = ps_type.upper()
    if ps_type == 'NCP':
        ps_suff = '_%s19_%s_OTF%s' % (ps_type, tmp, ps_ext)
    elif ps_type == 'SP':
        ps_suff = '_00%s' % ps_ext
    elif ps_type == 'EXT':
        ps_suff = '_%s_%s_OTF%s' % (ps_type, tmp, ps_ext)
    else:
        exit('problem with ps type...')

    for e in ele_list:
        calc.cell.species_pot = (e, e + ps_suff)

    return calc


def castep_calc_temp(template, xc_f="PBE", e_cut=800, kpts=[4, 4, 1], directory='data', ps_type='ext', ps_ext='.usp'):
    atoms = read(template)
    calc = atoms.calc

    # include interface settings in .param file
    calc._export_settings = True
    calc._pedantic = True

    # reuse the same directory
    calc._seed = directory
    calc._label = directory
    calc._directory = directory
    calc._rename_existing_dir = False

    # calc._link_pspots = True
    # calc._copy_pspots = False
    # calc._build_missing_pspots = False

    # Param settings
    calc.param.reuse = True
    calc.param.xc_functional = xc_f
    calc.param.cut_off_energy = e_cut
    calc.param.finite_basis_corr = 0
    calc.param.num_dump_cycles = 0  # Prevent CASTEP from writing *wvfn* files
    calc.param.write_checkpoint = 'minimal'
    calc.param.max_scf_cycles = 100
    calc.param.opt_strategy = "speed"
    calc.param.mixing_scheme = "Broyden"  # Broyden Pulay
    calc.param.spin_polarised = False
    # calc.param.elec_method = "EDFT"

    calc.param.sedc_apply = True
    calc.param.sedc_scheme = 'MBD*'
    calc.param.perc_extra_bands = 50
    calc.param.elec_energy_tol = 1e-6

    ele_list = atoms.get_chemical_symbols()
    tmp = None
    ps_suff = None
    if xc_f == 'PBE' or xc_f == 'PBE0':
        tmp = 'PBE'
    elif xc_f == 'B3LYP' or xc_f == 'BLYP':
        tmp = 'BLYP'
    elif xc_f == 'LDA':
        tmp = 'LDA'
    else:
        exit('problem with xc and not supporting ps type...')
    ps_type = ps_type.upper()
    if ps_type == 'NCP':
        ps_suff = '_%s19_%s_OTF%s' % (ps_type, tmp, ps_ext)
    elif ps_type == 'SP':
        ps_suff = '_00%s' % ps_ext
    elif ps_type == 'EXT':
        ps_suff = '_%s_%s_OTF%s' % (ps_type, tmp, ps_ext)
    else:
        exit('problem with ps type...')

    for e in ele_list:
        calc.cell.species_pot = (e, e + ps_suff)

    return calc


def make_graphene_sheet(i1=4, i2=4, alat=2.45, vaccum=12.8):
    # Set up a graphene lattice with a i1xi2 supercell
    return Graphene(symbol='C', latticeconstant={'a': alat, 'c': vaccum}, size=(i1, i2, 1))


def make_h_adsorbate(input_atoms, index=11, height=1.5, pos=None):
    if pos is None:
        # Get the site to adsorb the H atom
        pos = input_atoms.get_positions()
        pos = pos[index]
    # Make a hydrogen atom
    h_atom = Atoms('H', positions=[(pos[0], pos[1], pos[2] + height)])
    # Add a hydrogen atom to the slab
    return input_atoms + h_atom


f_max = 0.01
f_debug = None
f_optimise = False
f_mlneb = False
f_ts_refine = False

os_name = os.name
if os_name == 'nt':
    # Windows
    f_debug = True
else:
    # Linux
    f_debug = False

# gra = ase.io.castep.read_castep("4x4.castep")
reactant = read("4x4-out.cell")
# write('graphene_intial.traj', gra)
if f_debug:
    view(reactant)

# Apply the calculator
if f_debug:
    calc = EMT()
else:
    calc = castep_calc(reactant)
    #calc = castep_calc_temp("4x4-out.cell")

reactant.set_calculator(copy.deepcopy(calc))
t0 = time.time()
enegy = reactant.get_potential_energy()
t1 = time.time()
print('Time for energy calculation: %f s' % (t1 - t0))
print('Energy: %f eV' % enegy)

# # Optimise the inital state
# qn = BFGS(gra_s1, trajectory='graphene_opti.traj', logfile='opti.log')
# qn.run(fmax=f_max)
