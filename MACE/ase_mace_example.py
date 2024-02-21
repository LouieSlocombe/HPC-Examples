from mace.calculators import mace_anicc, mace_mp, mace_off
from ase.optimize import BFGS
from ase.visualize import view
from ase import build

# calc = mace_anicc()  # CC model nut requres a GPU
calc = mace_off(model="small", default_dtype="float64",device="cpu", dispersion=False)


# make a water molecule
atoms = build.molecule('H2O')
atoms.calc = calc

# optimise the geometry
BFGS(atoms).run(fmax=0.1)
view(atoms)