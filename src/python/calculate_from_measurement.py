from configurational_entropy import calculateConfigurationalEntropy
from measurement import Measurement

def calculateFromMeasurement(filepath):
	m = Measurement()
	m.fromMeasurement(filepath)
	n1, n2, H_n_values, c = m.readFile()
	calculateConfigurationalEntropy(n1, n2, H_n_values, c)
