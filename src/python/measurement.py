import os
from datetime import datetime

class Measurement():
	def __init__(self, filepath, covalent_radii_cut_off, c, n1, n2):
		self.filepath = filepath
		self.filename = os.path.basename(filepath)
		self.covalent_radii_cut_off = covalent_radii_cut_off
		self.c = c
		self.n1 = n1
		self.n2 = n2
		self.f = None
		self.createFile();

	def __init__(self, filepath):
		print("entrou aqui " + filepath) #TODO

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
