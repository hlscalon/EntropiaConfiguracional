import sys
import boost_graph as bg
import matplotlib.pyplot as plt
import os

from datetime import datetime
from ase import Atom
from ase.io import read
from ase.data import covalent_radii
from math import log
from operator import itemgetter
from random import uniform
from numpy import asarray
from numpy.polynomial.polynomial import polyfit

def run(G, m, n, slab, c):
	(pbcX, pbcY, pbcZ) = slab.get_pbc()
	(dmin, dmax) = getMaxMinSlabArray(slab)
	graphs = bg.Graphs(m, n, G, 3, len(slab), pbcX, pbcY, pbcZ, dmin, dmax)
	graphs.add_positions(slab.get_positions(wrap=True))
	cell = slab.get_cell()
	graphs.init_search(0, cell[0][0], 0, cell[1][1], 0, cell[2][2])
	Hc_n, valid = graphs.calculate_configurational_entropy(c)
	if not valid:
		print("n: %d. H1(n) exceeds 1%% of H(n). Not a valid measurement." % n)

	return Hc_n, valid

def generateSubgraphs(G, graphs, m, n, slab, atomic_number, aboria_tree):
	(dmin, dmax) = getMaxMinSlab(slab)

	closest_neighbors = []
	for i in range(m):
		(x, y, z) = generateRandomPoint(dmin, dmax)
		n_closest_neighbors = aboria_tree.search_nearest_neighbors(x, y, z, n)
		closest_neighbors.append(n_closest_neighbors)

	graphs.generate_subgraphs(G, n, closest_neighbors)

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

def main():
	if len(sys.argv) < 7:
		print("1 parameter: xyz filename\n2 parameter: covalent_radii_cut_off\n3 parameter: c\n4 parameter: initial n\n5 parameter: final n\n6 parameter: calculate (Y or N)")
		return

	filename = sys.argv[1]
	covalent_radii_cut_off = float(sys.argv[2]) # 1.12
	c = float(sys.argv[3])
	n1 = int(sys.argv[4])
	n2 = int(sys.argv[5])
	calculate = sys.argv[6] # Y or N

	if n1 > n2:
		print("Final m cannot be smaller than initial m")
		return

	print("Starting script...")

	slab = read(filename)

	print("Slab %s read with success" % filename)

	G = generateGraphFromSlab(slab, covalent_radii_cut_off)
	total_nodes = G.get_total_nodes()
	if total_nodes == 0 or G.get_total_edges() == 0:
		print("No edges found in graph. Check covalent_radii_cut_off")
		return

	print("Graph created with success. Nodes found: %d" % total_nodes)

	date_now = datetime.now()
	ce_file = "generated_files/gen_" + str(date_now.day) + "_" + str(date_now.month) + "_" + str(date_now.year) + "_" + str(date_now.hour) + "_" + str(date_now.minute) + "_" + str(date_now.second) + ".ce"
	dirname = os.path.dirname(ce_file)
	if not os.path.exists(dirname):
		os.makedirs(dirname)

	f = open(ce_file, "w+")
	f.write("filename: " + filename + "; covalent: " + str(covalent_radii_cut_off) + "; c: " + str(c) + "; n1: " + str(n1) + "; n2: " + str(n2) + "\r\n")

	hcn_values = []
	xy_polyfit = []

	for n in range(n1, n2):
		m = n * n * total_nodes
		(hcn, valid) = run(G, m, n, slab, c)

		f.write("n: " + str(n) + "; m: " + str(m) + "; hcn: " + str(hcn) + "; valid: " + str(valid) + "\r\n")

		hcn_values.append((n, hcn))
		if valid:
			xy_polyfit.append((n, hcn))

	f.close()

	if calculate == "Y":
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

	print("Program ended correctly")

if __name__ == "__main__":
	main()
