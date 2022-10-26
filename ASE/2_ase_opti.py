# Example: structure optimization of hydrogen molecule
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from ase.build import molecule
from ase.visualize import view
from ase.io import read

atoms = molecule('H2O')
view(atoms)
input("Press Enter to continue...")

atoms.calc = EMT()

opt = BFGS(atoms, trajectory='H2O.traj', logfile='H2O.log')
opt.run(fmax=0.01)


atom_opti = read("H2O.traj@:")

view(atom_opti)