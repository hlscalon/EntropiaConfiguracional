from ase.lattice.cubic import FaceCenteredCubic
from ase.visualize import view
from ase.io import write

atoms = FaceCenteredCubic(directions=[[1,-1,0], [1,1,-2], [1,1,1]],
                           size=(5,5,5), symbol='Cu', pbc=(1,1,1))

view(atoms)
write("fcc.xyz", atoms)
