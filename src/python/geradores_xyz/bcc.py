from ase.lattice.cubic import BodyCenteredCubic
from ase.visualize import view
from ase.io import write

atoms = BodyCenteredCubic(directions=[[1,0,0], [0,1,0], [1,1,1]],
						  size=(2,2,3), symbol='Cr', pbc=(1,1,1))

view(atoms)
write("bcc.xyz", atoms)
