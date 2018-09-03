from configurational_entropy import calculateConfigurationalEntropy
from measurement import Measurement

def calculateFromMeasurement(filepath):
	# calculateConfigurationalEntropy(n1, n2, xy_polyfit, hcn_values)
	m = Measurement(filepath) #TODO
	return