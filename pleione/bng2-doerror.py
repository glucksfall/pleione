# -*- coding: utf-8 -*-

'''
Project "Genetic Algorithm for Rule-Based Models", Rodrigo Santibáñez, 2017 @ Dlab, FCV (rsantibanez@dlab.cl)
A Python implementation of Alberto Martin's Genetic Algorithm, 2016 @ Dlab, FCV (ajmm@dlab.cl)
To be used with BNG2. Please refer to other subprojects for other stochastic simulators support
Citation:
'''

__author__  = 'Rodrigo Santibáñez'
__license__ = 'gpl-3.0'
__software__ = 'bng2-v2.3.2'

import argparse, io
import numpy, pandas
from pleione.fitness import do

def argsparser():
	parser = argparse.ArgumentParser(description = 'Calculate fitness between data and simulations. Only for BNG2 simulations.', \
		epilog = 'error acronysms are ADA, ANPWSD, APWSD, CHISQ, MNSE, MWUT, NPWSD, PWSD, SDA, SSQ, TOST, WMWET\n' \
			'see https://pleione.readthedocs.io/en/latest/ObjectiveFunctions.html for more information',
		formatter_class = argparse.RawTextHelpFormatter)

	# required args
	parser.add_argument('--data' , metavar = 'path' , type = str, required = True , nargs = '+', help = 'data (BNG2 format)')
	parser.add_argument('--sims' , metavar = 'path' , type = str, required = True , nargs = '+', help = 'BNG2 simulations')
	parser.add_argument('--file' , metavar = 'path' , type = str, required = True , nargs = 1  , help = 'output filename')
	parser.add_argument('--error', metavar = 'str'  , type = str, required = True , nargs = '+', help = 'fitness function(s) to calculate')

	# optional args
	# table of critical values for the Mann-Whitney U-test
	parser.add_argument('--crit' , metavar = 'path' , type = str, required = False, default = None, \
		help = 'Mann-Whitney U-test critical values')
	# report the matrices of the statistic tests
	parser.add_argument('--report', metavar = 'True', type = str, required = False, default = None, \
		help = 'report the arrays of the statistical tests')
	# calculate all fitness functions regardless of the used for model ranking
	parser.add_argument('--do_all', metavar = 'True', type = str, required = False, default = None, \
		help = 'calculate all fitness functions regardless of the used for ranking')

	return parser.parse_args()

# read simulation files
def read_sims(files):
	sims = []
	for infile in files:
		with open(infile, 'r') as file:
			tmp = io.StringIO(file.read()[1:]) # remove the # at the beginning of data files
			tmp = pandas.read_csv(tmp, delim_whitespace = True, header = 0, engine = 'python')
			tmp = tmp.set_index('time', drop = False).rename_axis(None, axis = 0).drop('time', axis = 1)
			sims.append(tmp)

	return pandas.concat(sims, keys = range(len(sims))), len(sims)

# read the data files
def read_data(files):
	data = []
	for infile in files:
		with open(infile, 'r') as file:
			tmp = io.StringIO(file.read()[1:]) # remove the # at the beginning of simulations files
			tmp = pandas.read_csv(tmp, delim_whitespace = True, header = 0, engine = 'python')
			tmp = tmp.set_index('time', drop = False).rename_axis(None, axis = 0).drop('time', axis = 1)
			data.append(tmp)

	return pandas.concat(data, keys = range(len(data))), len(data)

if __name__ == '__main__':
	args = argsparser()

	# read sims files
	sims, len_sims = read_sims(args.sims)
	# read data files
	data, len_data = read_data(args.data)

	# Filter out unavailable experimental data from simulation files and filter out non simulated observables
	sims = sims.filter(items = list(data.columns))
	data = data.filter(items = list(sims.columns))

	# Calculate fitness
	error = {}
	do(args, sims, len_sims, data, len_data, error)

	# write report file
	with open(args.file[0], 'w') as outfile:
		for fitfunc, value in sorted(error.items()):
			outfile.write('{:s}\t{:s}\n'.format(fitfunc, value))
