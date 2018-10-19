import sys

from src.python.configurational_entropy import startMeasurement
from src.python.calculate_from_measurement import calculateFromMeasurement

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
		if len(sys.argv) < 8:
			print("<$1>: from measurement (Y or N)\n" +
				  "<$2>: xyz filepath\n" +
				  "<$3>: file type (V = Vink file or N = Normal) \n" +
				  "<$4>: covalent radii cut off\n" +
				  "<$5>: c\n" +
				  "<$6>: initial n\n" +
				  "<$7>: final n\n" +
				  "<$8>: calculate (Y or N)")
			return

		filepath = sys.argv[2]
		file_type = sys.argv[3]
		covalent_radii_cut_off = float(sys.argv[4]) # 1.12
		c = float(sys.argv[5])
		n1 = int(sys.argv[6])
		n2 = int(sys.argv[7])
		calculate = sys.argv[8] # Y or N

		startMeasurement(filepath, file_type, covalent_radii_cut_off, c, n1, n2, calculate)

if __name__ == "__main__":
	main()
