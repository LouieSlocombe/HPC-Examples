from ase.visualize import view
from ase import Atoms

d = 1.1
co_atoms = Atoms('CO', positions=[(0, 0, 0), (0, 0, d)])
view(co_atoms)
input("Press Enter to continue...")

d = 2.9
L = 10.0
wire = Atoms('Au',
             positions=[[0, L / 2, L / 2]],  # position of the atom
             cell=[d, L, L],  # define the cell
             pbc=[1, 0, 0])  # make it periodic in x-direction
view(wire)
input("Press Enter to continue...")

from ase.build import molecule

atoms = molecule('H2O')
view(atoms)

from ase.io import write, read

atoms.rotate('z', 90, rotate_cell=True)
write('image.png', atoms)
write('H20.xyz', atoms)

saved = read('H20.xyz')

input("Press Enter to continue...")

# To setup an Al(111) surface with a hydrogen atom adsorbed in an on-top position
from ase.build import fcc111, add_adsorbate

slab = fcc111('Al', size=(2, 2, 3))
add_adsorbate(slab, 'H', 1.5, 'ontop')
slab.center(vacuum=20.0, axis=2)
view(slab)
