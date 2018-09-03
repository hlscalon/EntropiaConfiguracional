import os
from datetime import datetime

class Measurement():
	def __init__(self):
		return None

	def fromFile(self, filepath, covalent_radii_cut_off, c, n1, n2):
		self.filepath = filepath
		self.filename = os.path.basename(filepath)
		self.covalent_radii_cut_off = covalent_radii_cut_off
		self.c = c
		self.n1 = n1
		self.n2 = n2
		self.f = None

	def fromMeasurement(self, filepath):
		self.filepath = filepath

	def createFile(self):
		date_now = datetime.now()
		ce_file = "medicoes/med_" + self.filename + "_" + str(date_now.day) + "_" + str(date_now.month) + "_" + str(date_now.year) + "_" + str(date_now.hour) + "_" + str(date_now.minute) + "_" + str(date_now.second) + ".ce"
		dirname = os.path.dirname(ce_file)
		if not os.path.exists(dirname):
			os.makedirs(dirname)

		self.f = open(ce_file, "w+")
		self.f.write("filepath: " + self.filepath + "; covalent: " + str(self.covalent_radii_cut_off) + "; c: " + str(self.c) + "; n1: " + str(self.n1) + "; n2: " + str(self.n2) + "\r\n")

	def writeResult(self, n, m, hcn, valid):
		self.f.write("n: " + str(n) + "; m: " + str(m) + "; hcn: " + str(hcn) + "; valid: " + str(valid) + "\r\n")

	def close(self):
		self.f.close()

	def readFile(self):
		f = open(self.filepath, "r")

		firstLine = f.readline()
		n1idx = firstLine.find("n1: ")
		n2idx = firstLine.find("; n2: ")

		n1 = int(firstLine[n1idx + 4 : n2idx])
		n2 = int(firstLine[n2idx + 6 : ])

		xy_polyfit = []
		hcn_values = []
		for line in f:
			nidx = line.find("n: ")
			midx = line.find("; m: ")
			n = int(line[nidx + 3 : midx])
			hcnidx = line.find("hcn: ")
			valididx = line.find("; valid: ")
			hcn = float(line[hcnidx + 5 : valididx])
			valid = line[valididx + 9 : ].strip() == 'True'

			hcn_values.append((n, hcn))
			if valid:
				xy_polyfit.append((n, hcn))

		f.close()

		return n1, n2, xy_polyfit, hcn_values
