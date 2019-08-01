# -*- coding: utf-8 -*-

'''
Project "Genetic Algorithm for Rule-Based Models", Rodrigo Santibáñez, 2017 @ Dlab, FCV (rsantibanez@dlab.cl)
A Python implementation of Alberto Martin's Genetic Algorithm, 2016 @ Dlab, FCV (ajmm@dlab.cl)
To be used with PISKaS. Please refer to other subprojects for other stochastic simulators support
Citation:
'''

__author__  = 'Rodrigo Santibáñez'
__license__ = 'gpl-3.0'
__software__ = 'piskas-v1.3'

import argparse, sys
import pandas, numpy

def argsparser():
	parser = argparse.ArgumentParser(description = 'Calculate goodness of fit between data and simulations.')
	parser.add_argument('--data' , metavar = 'path', type = str, required = True , nargs = '+', help = 'data files')
	parser.add_argument('--sims' , metavar = 'path', type = str, required = True , nargs = '+', help = 'PISKaS output without further processing')
	parser.add_argument('--file' , metavar = 'path', type = str, required = True , nargs = 1  , help = 'output file name')
	parser.add_argument('--crit' , metavar = 'path', type = str, required = False, nargs = 1  , help = 'Mann-Whitney U-test critical values')

	return parser.parse_args()

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
	args = argsparser()

	# read sims files
	sims, len_sims = read_sims(args.sims)
	# read data files
	data, len_data = read_data(args.data)

	# Filter out unavailable experimental data from simulation files and filter out unsimulated observables
	sims = sims.filter(items = list(data.columns))
	data = data.filter(items = list(sims.columns))

	# Calculate fitness
	error = {}
	do(args, sims, len_sims, data, len_data, error, True)

	# write report file
	with open(args.file[0], 'w') as outfile:
		for fitfunc, value in sorted(error.items()):
			outfile.write('{:s}\t{:s}\n'.format(fitfunc, value))
