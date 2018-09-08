import sys
import boost_graph as bg
import matplotlib.pyplot as plt

from measurement import Measurement
from ase import Atom
from ase.io import read
from ase.data import covalent_radii
from ase.geometry import is_orthorhombic
from math import log
from operator import itemgetter
from random import uniform
from numpy import asarray
from numpy.polynomial.polynomial import polyfit

def run(configurationEntropy, m, n, c):
	Hc_n, valid = configurationEntropy.calculate(m, n, c)
	if not valid:
		print("n: %d. H1(n) exceeds 1%% of H(n). Not a valid measurement." % n)

	return Hc_n, valid

def getMaxMinSlabArray(slab):
	(dmin, dmax) = (
		[float("Inf"), float("Inf"), float("Inf")],  # x, y, z
		[-float("Inf"), -float("Inf"), -float("Inf")]   # x, y, z
	)

	positions = slab.get_positions(wrap=True) # wrap atoms back to simulation cell
	for distance in positions:
		for idx, d in enumerate(distance):
			if (d > dmax[idx]):
				dmax[idx] = d
			if (d < dmin[idx]):
				dmin[idx] = d

	return (asarray(dmin), asarray(dmax))

def generateGraphFromSlab(slab, covalent_radii_cut_off):
	graph = bg.Graph()

	atomic_numbers = slab.get_atomic_numbers()
	all_distances = slab.get_all_distances(mic=True)
	for atom1, distances in enumerate(all_distances):
		if not graph.has_node(atom1):
			graph.add_node(atom1) # add nodes not bonded

		atom1_cr = covalent_radii[atomic_numbers[atom1]]
		for atom2, distance in enumerate(distances):
			if atom1 != atom2:
				atom2_cr = covalent_radii[atomic_numbers[atom2]]
				# if the distance between two atoms is less than the sum of their covalent radii, they are considered bonded.
				if (distance < ((atom1_cr + atom2_cr) * covalent_radii_cut_off)):
					graph.add_edge(atom1, atom2)

	return graph

def calculateConfigurationalEntropy(n1, n2, xy_polyfit, hcn_values):
	(x_p, y_p) = zip(*xy_polyfit)
	x_p = asarray(x_p)
	y_p = asarray(y_p)

	# straight line fit
	b, m = polyfit(x_p, y_p, 1) # m equals the slope of the line
	plt.plot(x_p, b + m * x_p, '-')

	x, y = zip(*hcn_values)
	plt.scatter(x, y)

	plt.axis([n1, n2, -5, 10])
	plt.show()

	print("Estimated configurational entropy = %f" % (m))


def startMeasurement(filepath, covalent_radii_cut_off, c, n1, n2, calculate):
	if n1 > n2:
		print("Final m cannot be smaller than initial m")
		return

	print("Starting script...")

	slab = read(filepath)

	print("Slab %s read with success" % filepath)

	G = generateGraphFromSlab(slab, covalent_radii_cut_off)
	total_nodes = G.get_total_nodes()
	if total_nodes == 0 or G.get_total_edges() == 0:
		print("No edges found in graph. Check covalent_radii_cut_off")
		return

	print("Graph created with success. Nodes found: %d" % total_nodes)

	measurement = Measurement()
	measurement.fromFile(filepath, covalent_radii_cut_off, c, n1, n2)
	measurement.createFile()

	hcn_values = []
	xy_polyfit = []

	(pbcX, pbcY, pbcZ) = slab.get_pbc()
	(dmin, dmax) = getMaxMinSlabArray(slab)
	configurationalEntropy = bg.ConfigurationalEntropy(G, 3, len(slab), pbcX, pbcY, pbcZ, dmin, dmax)
	configurationalEntropy.add_positions(slab.get_positions(wrap=True))
	cell = slab.get_cell()

	if not is_orthorhombic(cell):
		print("Unit cell is not orthorhombic")
		return

	# somente para celulas ortogonais
	configurationalEntropy.init_search(0, cell[0][0], 0, cell[1][1], 0, cell[2][2])

	for n in range(n1, n2):
		m = n * n * total_nodes
		(hcn, valid) = run(configurationalEntropy, m, n, c)

		measurement.writeResult(n, m, hcn, valid)

		hcn_values.append((n, hcn))
		if valid:
			xy_polyfit.append((n, hcn))

	measurement.close()

	if calculate == "Y":
		calculateConfigurationalEntropy(n1, n2, xy_polyfit, hcn_values)

	print("Program ended correctly")
