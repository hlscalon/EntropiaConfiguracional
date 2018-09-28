import sys

from src.python.configurational_entropy import startMeasurement
from src.python.calculate_from_measurement import calculateFromMeasurement

"""
python2 main.py N arquivos_xyz/SILICA_SAMPLES/3k/3k-si-si.xyz 1.62 0 3 10 Y
python2 main.py N arquivos_xyz/SILICA_SAMPLES/3k/si-si-0-17-ase.xyz 1.62 0 3 20 Y
python2 main.py N arquivos_xyz/SILICA_SAMPLES/3k/3k-0-17-ase.xyz 1.12 0 3 15 Y
python2 main.py N arquivos_xyz/SILICA_SAMPLES/3k/3k-0-12-ase.xyz 1.12 0 3 15 Y
python2 main.py N arquivos_xyz/fcc.xyz 1.12 0 3 10 Y
python2 main.py N arquivos_xyz/fcc.xyz 1.12 0 3 15 Y
time python2 main.py N arquivos_xyz/fcc.xyz 1.12 0 3 30 Y
python2 main.py Y medicoes/med_fcc.xyz_3_9_2018_7_58_35.ce
python2 -m cProfile -s time main.py N arquivos_xyz/fcc.xyz 1.12 0 3 10 Y
python2 -m yep --  main.py N arquivos_xyz/fcc.xyz 1.12 0 3 15 Y
pprof --web /usr/bin/python2.7-debug main.py.prof
"""

def main():
	if len(sys.argv) < 2:
		print("1 parameter: from measurement (Y or N)\n")
		return

	measurement = sys.argv[1]

	if measurement == "Y":
		if len(sys.argv) < 3:
			print("1 parameter: from measurement (Y or N)\n" +
				  "2 parameter: measurement filepath")
			return

		filepath = sys.argv[2]

		calculateFromMeasurement(filepath)
	else:
		if len(sys.argv) < 7:
			print("1 parameter: from measurement (Y or N)\n" +
				  "2 parameter: xyz filepath\n" +
				  "3 parameter: covalent_radii_cut_off\n" +
				  "4 parameter: c\n" +
				  "5 parameter: initial n\n" +
				  "6 parameter: final n\n" +
				  "7 parameter: calculate (Y or N)")
			return

		filepath = sys.argv[2]
		covalent_radii_cut_off = float(sys.argv[3]) # 1.12
		c = float(sys.argv[4])
		n1 = int(sys.argv[5])
		n2 = int(sys.argv[6])
		calculate = sys.argv[7] # Y or N

		startMeasurement(filepath, covalent_radii_cut_off, c, n1, n2, calculate)

if __name__ == "__main__":
	main()
