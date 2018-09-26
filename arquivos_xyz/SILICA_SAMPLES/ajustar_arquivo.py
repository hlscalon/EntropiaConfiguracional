import ase
import networkx as nx
import matplotlib.pyplot as plt

from ase.io import read, write
from ase.visualize import view
from ase.data import covalent_radii

slab = read("3k/3k-0-12-ase.xyz")
view(slab)
# slab.wrap()
# view(slab)
# slab.write("3k/3k-ase.xyz")

"""
Si-O = rc 1.8
Si-O-Si = rc 3.6
Si-Si = rc 1.11

cutoff = 3.6 -> 2.22 = 1.62
"""

# criar arquivo com atomos dentro da cell
"""
cell = slab.get_cell()
novo_slab = ase.Atoms() #slab
novo_slab.set_cell(cell)
novo_slab.set_pbc(True)
for idx, atom in enumerate(slab):
	if atom.x >= 0 and atom.x <= cell[0][0] and \
	   atom.y >= 0 and atom.y <= cell[1][1] and \
	   atom.z >= 0 and atom.z <= cell[2][2]:
		novo_slab.append(atom)

view(novo_slab)
novo_slab.write("3k/3k-0-10-ase.xyz")
"""

# verificar grafo gerado
#"""
covalent_radii_cut_off = 1.12
all_distances = slab.get_all_distances(mic=True)
atomic_numbers = slab.get_atomic_numbers()

graph = nx.Graph()
"""
for atom1, distances in enumerate(all_distances):
	if not atom1 in graph:
		graph.add_node(atom1) # add nodes not bonded

	atom1_cr = covalent_radii[atomic_numbers[atom1]]
	for atom2, distance in enumerate(distances):
		if atom1 != atom2:
			atom2_cr = covalent_radii[atomic_numbers[atom2]]
			# if the distance between two atoms is less than the sum of their covalent radii, they are considered bonded.
			if (distance < ((atom1_cr + atom2_cr) * covalent_radii_cut_off)):
				graph.add_edge(atom1, atom2)
"""
mapping = {}
for atom1, distances in enumerate(all_distances):
	if atomic_numbers[atom1] == 14: # Si
		atom1_cr = covalent_radii[atomic_numbers[atom1]]
		for atom2, distance in enumerate(distances):
			if atom1 != atom2:
				if atomic_numbers[atom2] == 8:
					atom2_cr = covalent_radii[atomic_numbers[atom2]]
					# if the distance between two atoms is less than the sum of their covalent radii, they are considered bonded.
					if (distance < ((atom1_cr + atom2_cr) * covalent_radii_cut_off)):
						try:
							mapping[atom1][atom2] = True
						except KeyError:
							mapping[atom1] = {}
							mapping[atom1][atom2] = True
							pass

						try:
							mapping[atom2][atom1] = True
						except KeyError:
							mapping[atom2] = {}
							mapping[atom2][atom1] = True
							pass

for atom1, at1 in enumerate(slab):
	if atomic_numbers[atom1] == 14: # Si
		if not atom1 in graph:
			graph.add_node(atom1) # add nodes not bonded

		for n1_o in mapping[atom1]:
			if atomic_numbers[n1_o] == 8: # O
				for n2_si in mapping[n1_o]:
					if atomic_numbers[n2_si] == 14: # Si
						graph.add_edge(atom1, n2_si)

#
nx.draw(graph, with_labels=True, font_weigth='bold')
plt.show()

# remove O atomos
del slab[[atom.index for atom in slab if atom.symbol == 'O']]

#slab.write("3k/3k-0-12-ase-2.xyz")
print slab
#"""
