# -*- coding: utf-8 -*-

'''
Project "Genetic Algorithm for Rule-Based Models", Rodrigo Santibáñez, 2017 @ Dlab, FCV (rsantibanez@dlab.cl)
A Python implementation of Alberto Martin's Genetic Algorithm, 2016 @ Dlab, FCV (ajmm@dlab.cl)
To be used with KaSim v3.5. Please refer to other subprojects for other stochastic simulators support
Citation:
'''

__author__  = 'Rodrigo Santibáñez'
__license__ = 'gpl-3.0'
__software__ = 'kasim-v3.5'

import argparse, os
import pandas, numpy
from fitness import do

def argsparser():
	parser = argparse.ArgumentParser(description = 'Calculate goodness of fit between data and simulations.')
	parser.add_argument('--data' , metavar = 'path', type = str, required = True , nargs = '+', help = 'data files')
	parser.add_argument('--sims' , metavar = 'path', type = str, required = True , nargs = '+', help = 'KaSim v3.5 output without further processing')
	parser.add_argument('--file' , metavar = 'path', type = str, required = True , nargs = 1  , help = 'output file name')
	parser.add_argument('--error', metavar = 'str' , type = str, required = True , nargs = '+', help = 'Goodness of Fit Function(s) to calculate')
	parser.add_argument('--crit' , metavar = 'path', type = str, required = False, nargs = 1  , help = 'Mann-Whitney U-test critical values')

	# DEPRECATED path to R executable and libs
	#parser.add_argument('--r_path', metavar = 'path', type = str, required = False, default = '~/bin/R', help = 'R exe path, default ~/bin/R')
	#parser.add_argument('--r_libs', metavar = 'path', type = str, required = False, default = ''       , help = 'R lib path, default empty')
	# report MWUT, WMWET?
	parser.add_argument('--report', metavar = 'str' , type = str, required = False, default = None     , help = 'report the array of U-tests and/or Wellek\'s tests')

	return parser.parse_args()

# read simulation files
def read_sims(files):
	sims = []
	for infile in files:
		with open(infile, 'r') as file:
			sims.append(pandas.read_csv(file, delimiter = ',', skiprows = 2, header = 0, engine = 'python').set_index('[T]', drop = False).rename_axis(None, axis = 0).drop('[T]', axis = 1))

	return pandas.concat(sims, keys = range(len(sims))), len(sims)

# read the data files
def read_data(files):
	data = []
	for infile in files:
		with open(infile, 'r') as file:
			data.append(pandas.read_csv(file, delimiter = ',', header = 0, engine = 'python').set_index('[T]', drop = False).rename_axis(None, axis = 0).drop('[T]', axis = 1))

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
	do(args, sims, len_sims, data, len_data, error, False)

	# write report file
	with open(args.file[0], 'w') as outfile:
		for fitfunc, value in sorted(error.items()):
			outfile.write('{:s}\t{:s}\n'.format(fitfunc, value))
