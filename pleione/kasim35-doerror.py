# -*- coding: utf-8 -*-

'''
Project "Genetic Algorithm for Rule-Based Models", Rodrigo Santibáñez, 2017 @ Dlab, FCV (rsantibanez@dlab.cl)
A Python implementation of Alberto Martin's Genetic Algorithm, 2016 @ Dlab, FCV (ajmm@dlab.cl)
To be used with KaSim v3.5. Please refer to other subprojects for other stochastic simulators support
Citation: Pleione: A tool for statistical and multi-objective calibration of Rule-based models. Scientific Reports (2019)
DOI:
'''

__author__  = 'Rodrigo Santibáñez'
__license__ = 'gpl-3.0'
__software__ = 'kasim-v3.5'

import numpy, pandas
from pleione.fitness import argsparser, doerror

# read simulation files
def read_sims(files):
	sims = []
	for infile in files:
		with open(infile, 'r') as file:
			tmp = pandas.read_csv(file, delimiter = ' ', header = 0, engine = 'python')
			tmp = tmp.set_index('time', drop = False).rename_axis(None, axis = 0).drop('time', axis = 1)
			sims.append(tmp)

	return pandas.concat(sims, keys = range(len(sims))), len(sims)

# read the data files
def read_data(files):
	data = []
	for infile in files:
		with open(infile, 'r') as file:
			tmp = pandas.read_csv(file, delimiter = ' ', header = 0, engine = 'python')
			tmp = tmp.set_index('time', drop = False).rename_axis(None, axis = 0).drop('time', axis = 1)
			data.append(tmp)

	return pandas.concat(data, keys = range(len(data))), len(data)

if __name__ == '__main__':
	args = argsparser(**{ 'simulator' : 'KaSim v3.5'})
	data, len_data = read_data(args.data) # read data files
	sims, len_sims = read_sims(args.sims) # read sims files
	# calculate fitness
	doerror(args, data, len_data, sims, len_sims)
