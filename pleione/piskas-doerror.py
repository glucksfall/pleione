# -*- coding: utf-8 -*-

'''
Project "Genetic Algorithm for Rule-Based Models", Rodrigo Santibáñez, 2017 @ Dlab, FCV (rsantibanez@dlab.cl)
A Python implementation of Alberto Martin's Genetic Algorithm, 2016 @ Dlab, FCV (ajmm@dlab.cl)
To be used with PISKaS. Please refer to other subprojects for other stochastic simulators support
Citation: Pleione: A tool for statistical and multi-objective calibration of Rule-based models. Scientific Reports (2019)
DOI:
'''

__author__  = 'Rodrigo Santibáñez'
__license__ = 'gpl-3.0'
__software__ = 'piskas-v1.3'

import numpy, pandas
from pleione.fitness import argsparser, doerror

# read simulation files
def read_sims(files):
	sims = []
	for ind, infile in enumerate(files):
		with open(infile, 'r') as input:
			sims.append(pandas.read_csv(input, delimiter = ' ', engine = 'python', skipfooter = 9).set_index('time', drop = True))

		# mark each sim file with the corresponding compartment
		name = infile.split('.')[-2]
		sims[ind]['compartment'] = name
		sims[ind].set_index('compartment', append = True, inplace = True)
		sims[ind] = sims[ind].reorder_levels(['compartment', 'time'])

	return pandas.concat(sims, keys = range(len(sims))), len(sims)

# read the data files
def read_data(files):
	data = []
	for ind, infile in enumerate(files):
		with open(infile, 'r') as input:
			data.append(pandas.read_csv(input, delimiter = ' ', engine = 'python').set_index('time'))

		# mark each data file with the corresponding compartment
		name = list(data[ind].columns)[0]
		data[ind]['compartment'] = name
		data[ind].set_index('compartment', append = True, inplace = True)
		data[ind] = data[ind].reorder_levels(['compartment', 'time'])

	return pandas.concat(data, keys = range(len(data))), len(data)

if __name__ == '__main__':
	args = argsparser(**{ 'simulator' : 'PISKaS v1.3'})
	data, len_data = read_data(args.data) # read data files
	sims, len_sims = read_sims(args.sims) # read sims files
	# calculate fitness
	doerror(args, data, len_data, sims, len_sims)
