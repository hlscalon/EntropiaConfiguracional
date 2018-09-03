from configurational_entropy import calculateConfigurationalEntropy
from measurement import Measurement

def calculateFromMeasurement(filepath):
	m = Measurement()
	m.fromMeasurement(filepath)
	n1, n2, xy_polyfit, hcn_values = m.readFile()
	calculateConfigurationalEntropy(n1, n2, xy_polyfit, hcn_values)
