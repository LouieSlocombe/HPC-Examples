# Example: structure optimization of hydrogen molecule
from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.nwchem import NWChem
from ase.io import write
from ase.build import molecule
from ase.visualize import view
from ase.io import write, read

atoms = molecule('H2O')

calc = NWChem(label='calc/nwchem',
              dft=dict(maxiter=100,
                       xc='B3LYP'),
              basis='6-31+G*')

atoms.calc = calc

opt = BFGS(atoms, trajectory='H2O.traj', logfile='H2O.log')
opt.run(fmax=0.01)

atom_opti = read("H2O.traj@:")
