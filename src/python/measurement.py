import os
from datetime import datetime
from timeit import default_timer

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
		self.ce_file = ""
		self.iniTime = 0
		self.endTime = 0

	def fromMeasurement(self, filepath):
		self.filepath = filepath

	def createFile(self):
		date_now = datetime.now()
		self.ce_file = "medicoes/med_" + self.filename + "_" + str(date_now.day) + "_" + str(date_now.month) + "_" + str(date_now.year) + "_" + str(date_now.hour) + "_" + str(date_now.minute) + "_" + str(date_now.second) + ".ce"
		dirname = os.path.dirname(self.ce_file)
		if not os.path.exists(dirname):
			os.makedirs(dirname)

		with open(self.ce_file, "w+") as file:
			file.write("filepath: " + self.filepath + "; covalent: " + str(self.covalent_radii_cut_off) + "; c: " + str(self.c) + "; n1: " + str(self.n1) + "; n2: " + str(self.n2) + ";\n")

	def writeResult(self, n, m, H_n, H1n):
		line = "n: " + str(n) + "; m: " + str(m) + "; H_n: " + str(H_n) + "; H1n: " + str(H1n) + "; consider: Y;\n"

		with open(self.ce_file, "a") as file:
			file.write(line)

		print line,

	def start(self):
		self.iniTime = default_timer()

	def end(self):
		self.endTime = default_timer()
		print "Time taken: " + str(self.endTime - self.iniTime)

	def readFile(self):
		f = open(self.filepath, "r")

		firstLine = f.readline()
		cidx = firstLine.find("; c: ")
		n1idx = firstLine.find("; n1: ")
		n2idx = firstLine.find("; n2: ")
		fimidx = firstLine.find(";\n")

		c = float(firstLine[cidx + 5 : n1idx])
		n1 = int(firstLine[n1idx + 6 : n2idx])
		n2 = int(firstLine[n2idx + 6 : fimidx])

		H_n_values = []
		for line in f:
			nidx = line.find("n: ")
			midx = line.find("; m: ")
			n = int(line[nidx + 3 : midx])
			H_nidx = line.find("H_n: ")
			H1nidx = line.find("; H1n: ")
			consideridx = line.find("; consider: ")
			fimidx = line.find(";\n")
			H_n = float(line[H_nidx + 5 : H1nidx])
			H1n = float(line[H1nidx + 7 : consideridx])
			consider = line[consideridx + 12 : fimidx]

			H_n_values.append((n, H_n, H1n, consider))

		f.close()

		return n1, n2, H_n_values, c
