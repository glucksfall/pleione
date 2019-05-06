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

def do(error):
	"""
	# Fitness Calculation Template
	func = 0
	func = an algebraic expression combining the data average (data_avrg), data variance (data_stdv), simulation average (sims_stdv),
	single experimental files (data.loc[i]) and/or simulation files (sims.loc[i]).
	# Please consider this variables are DataFrames, meaning that division is a method (pandas.DataFrame.division)
	# drop NaN values (from experimental data without simulation point or vice-versa), sum the two dimensions, and return a 6 float points scientific notation number
	error['acronysm'] = '{:.6e}'.format(func.dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all').sum().sum())
	"""

	# Calculate common variables this script use further to measure fitness
	data_avrg = 0
	for i in range(len_data):
		data_avrg += data.loc[i].divide(len_data)

	data_stdv = 0
	for i in range(len_data):
		data_stdv += ((data.loc[i] - data_avrg)**2).divide(len_data - 1)
	data_stdv = data_stdv.filter(items = list(sims.columns))

	sims_avrg = 0
	for i in range(len_sims):
		sims_avrg += sims.loc[i].divide(len_sims)

	# mean square error
	func = 0
	func = (data_avrg - sims_avrg)**2
	error['SDM'] = '{:.6e}'.format(func.dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all').sum().sum())

	# mean absolute error
	func = 0
	func = abs(data_avrg - sims_avrg)
	error['MAE'] = '{:.6e}'.format(func.dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all').sum().sum())

	# sum of squares (from BioNetFit paper)
	func = 0
	for i in range(len_data):
		for j in range(len_sims):
			func += (data.loc[i] - sims.loc[j])**2
	error['SSQ'] = '{:.6e}'.format(func.dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all').sum().sum())

	# mean normalized sum of squares (from BioNetFit paper)
	func = 0
	for i in range(len_data):
		for j in range(len_sims):
			func += ((data.loc[i] - sims.loc[j]).divide(data_avrg))**2
	error['MNSE'] = '{:.6e}'.format(func.replace([numpy.inf, -numpy.inf], numpy.nan).dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all').sum().sum())

	# chi-square (from BioNetFit paper)
	func = 0
	for i in range(len_data):
		for j in range(len_sims):
			func += ((data.loc[i] - sims.loc[j]).divide(data_stdv**0.5))**2
	error['CHISQ'] = '{:.6e}'.format(func.dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all').sum().sum())

	# pair-wise square deviation
	func = 0
	for i in range(len_data):
		for j in range(len_sims):
			func += ((data.loc[i] - sims.loc[j])**2).divide(len_data * len_sims)
	error['PWSD'] = '{:.6e}'.format(func.dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all').sum().sum())

	# pair-wise absolute deviation
	func = 0
	for i in range(len_data):
		for j in range(len_sims):
			func += (abs(data.loc[i] - sims.loc[j])).divide(len_data * len_sims)
	error['APWSD'] = '{:.6e}'.format(func.dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all').sum().sum())

	# normalized pair-wise square deviation (also implemented in BioNetFit as equation 3, but not normalized by the number of data * sims)
	func = 0
	for i in range(len_data):
		for j in range(len_sims):
			func += (((data.loc[i] - sims.loc[j]).divide(data.loc[i]))**2).divide(len_data * len_sims)
	error['NPWSD'] = '{:.6e}'.format(func.replace([numpy.inf, -numpy.inf], numpy.nan).dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all').sum().sum())

	# normalized pair-wise absolute deviation
	func = 0
	for i in range(len_data):
		for j in range(len_sims):
			func += (abs((data.loc[i] - sims.loc[j]).divide(data.loc[i]))).divide(len_data * len_sims)
	error['ANPWSD'] = '{:.6e}'.format(func.replace([numpy.inf, -numpy.inf], numpy.nan).dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all').sum().sum())

	# Mann-Whitney U-test
	if ((len_data >= 3 and len_sims >= 3) or (len_data >= 2 and len_sims >= 5)):

		ucrit = pandas.read_table(args.crit[0], sep = None, engine = 'python', header = 0, index_col = 0)
		udata = pandas.DataFrame(index = sims.loc[0].index, columns = sims.loc[0].columns).fillna(0)
		usims = pandas.DataFrame(index = sims.loc[0].index, columns = sims.loc[0].columns).fillna(0)

		for i in range(len_data):
			for j in range(len_sims):
				diff = (data.loc[i] - sims.loc[j]).dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all')
				# if data < sims count -1.0
				diff[diff < 0] = -1.0
				# if data > sims count +1.0
				diff[diff > 0] = +1.0
				# if data = sims count +0.5
				diff[diff == 0] = +0.5
				# count how many times is data < sims (udata and usims are complementary)
				udata += diff[diff == -1.0].fillna(0).divide(-1) + diff[diff == +0.5].fillna(0)
				usims += diff[diff == +1.0].fillna(0).divide(+1) + diff[diff == +0.5].fillna(0)

		# U is significant if it is less than or equal to the table value
		U = len_data * len_sims - udata.where(udata >= usims).fillna(usims.where(usims >= udata))
		U[U <= ucrit.loc[len_sims, str(len_data)]] = +1.0
		U[U > ucrit.loc[len_sims, str(len_data)]] = +0.0

		error['MWUT'] = '{:.0f}'.format(U.sum().sum())
	else:
		error['MWUT'] = str(numpy.nan)

if __name__ == '__main__':
	args = argsparser()

	# read sims files
	sims, len_sims = read_sims(args.sims)
	# read data files
	data, len_data = read_data(args.data)

	# Filter out unavailable experimental data from simulation files and filter out unsimulated observables
	sims = sims.filter(items = list(data.columns))
	data = data.filter(items = list(sims.columns))

	# add here your favorite error function and call the genetic algorithm script with its acronysm
	error = {
		'SDM'   : str(numpy.nan),
		'MAE'   : str(numpy.nan),
		'SSQ'   : str(numpy.nan),
		'MWUT'  : str(numpy.nan),
		'PWSD'  : str(numpy.nan),
		'MNSE'  : str(numpy.nan),
		'CHISQ' : str(numpy.nan),
		'APWSD' : str(numpy.nan),
		'NPWSD' : str(numpy.nan),
		'ANPWSD': str(numpy.nan),
		}

	do(error)

	# write report file
	with open(args.file[0], 'w') as outfile:
		for fitfunc, value in sorted(error.items()):
			outfile.write('{:s}\t{:s}\n'.format(fitfunc, value))
