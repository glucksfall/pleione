# -*- coding: utf-8 -*-

'''
Project "Genetic Algorithm for Rule-Based Models", Rodrigo Santibáñez, 2017 @ Dlab, FCV (rsantibanez@dlab.cl)
A Python implementation of Alberto Martin's Genetic Algorithm, 2016 @ Dlab, FCV (ajmm@dlab.cl)
Citation: Pleione: A tool for statistical and multi-objective calibration of Rule-based models. Scientific Reports (2019)
DOI:
'''

__author__  = 'Rodrigo Santibáñez'
__license__ = 'gpl-3.0'

import argparse
import numpy, pandas

# argparser
def argsparser(**kwargs):
	parser = argparse.ArgumentParser(description = 'Calculate fitness between data and simulations. Only for {:s} simulations'.format(kwargs['simulator']), \
		epilog = 'error acronysms are ADA, ANPWSD, APWSD, CHISQ, DUT, MNSE, MWUT, NPWSD, PWSD, SDA, SSQ, TOST, WMWET\n' \
			'see https://pleione.readthedocs.io/en/latest/ObjectiveFunctions.html for more information',
		formatter_class = argparse.RawTextHelpFormatter)

	# required args
	parser.add_argument('--data'  , metavar = 'path' , type = str, required = True , nargs = '+', help = 'data ({:s} format)'.format(kwargs['simulator']))
	parser.add_argument('--sims'  , metavar = 'path' , type = str, required = True , nargs = '+', help = '{:s} simulations'.format(kwargs['simulator']))
	parser.add_argument('--file'  , metavar = 'path' , type = str, required = True , nargs = 1  , help = 'output filename')
	parser.add_argument('--error' , metavar = 'str'  , type = str, required = False, nargs = '+', help = 'fitness function(s) to calculate')

	#subparsers = parser.add_subparsers()
	# optional args
	#opts = subparsers.add_parser(name = 'a', help = 'report all calculations')
	# report the matrices of the statistic tests
	parser.add_argument('--report', metavar = 'True' , type = str, required = False, default = None, \
		help = 'report the arrays of the statistical tests')
	# calculate all fitness functions regardless of the used for model ranking
	parser.add_argument('--do_all', metavar = 'True' , type = str, required = False, default = None, \
		help = 'calculate all fitness functions regardless of the used for ranking')

	# more optional args (for equivalence tests)
	#equiv = subparsers.add_parser(name = 'b', help = 'optional for equivalence tests')
	parser.add_argument('--crit'  , metavar = 'path' , type = str, required = False, default = None, \
		help = 'Mann-Whitney U-test critical values')
	parser.add_argument('--lower' , metavar = 'path' , type = str, required = False, default = None, \
		help = 'file with the lower limit for the equivalence test. Same format as data')
	parser.add_argument('--upper' , metavar = 'path' , type = str, required = False, default = None, \
		help = 'file with the upper limit for the equivalence test. Same format as data\n' \
			'Setting either lower or upper will make the threshold symmetric.')
	parser.add_argument('--stdv'  , metavar = 'sims' , type = str, required = False, default = 'data', \
		help = 'use the simulation standard deviation (sims stdv) instead of data stdv as lower and upper as limits')
	parser.add_argument('--factor', metavar = 'float', type = str, required = False, default = '1' , \
		help = 'factor to divide lower and upper in case of using data stdv or sims stdv.')

	return parser.parse_args()

# main
def doerror(args, data, len_data, sims, len_sims):
	# Filter out unavailable experimental data from simulation files and filter out non simulated observables
	sims = sims.filter(items = list(data.columns))
	data = data.filter(items = list(sims.columns))

	# Calculate fitness
	error = {}
	docalc(args, data, len_data, sims, len_sims, error)

	# write report file
	with open(args.file[0], 'w') as outfile:
		for fitfunc, value in sorted(error.items()):
			outfile.write('{:s}\t{:s}\n'.format(fitfunc, value))

	return 0

# helpers
def doavrg(data, len_data):
	avrg = 0
	for i in range(len_data):
		avrg += data.loc[i].divide(len_data)
	return avrg

def dostdv(data, len_data):
	avrg = doavrg(data, len_data)
	stdv = 0
	if len_data > 1:
		for i in range(len_data):
			stdv += ((data.loc[i] - avrg)**2).divide(len_data - 1).replace(0., 1.) # replications with stdv = 0 are problematic, because stdv is used to divide a number
	else:
		stdv = pandas.DataFrame(index = data.loc[0].index, columns = data.loc[0].columns).fillna(1.) # a DataFrame which every stdv is equal to 1.
	return stdv**.5

# calc error
def docalc(args, data, len_data, sims, len_sims, error):
	"""
	# Fitness Calculation Template:
	if set(args.error).issuperset(set(['the-acronysm'])):
		1. func = 0

		2. func = an algebraic expression combining the data average (data_avrg), data standard deviation (data_stdv), simulation average (sims_stdv),
		simulation standard deviation (sims_stdv), single experimental files (data.loc[i]), and/or simulation files (sims.loc[j])
		Note1: Perform two for-loops if using data.loc[i] and sims.loc[j].
		Note2: Please consider these variables are DataFrames, meaning that multiplication and division are methods (e.g. df1.division(df2))

		3. Drop NaN values (from experimental time points without simulated values, or simulated values without experimental data)
		with dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all'). Also transform Inf values with replace([numpy.inf, -numpy.inf], numpy.nan)

		4. Sum the two dimensions, and return a 6 float points scientific notation number (0 float points for statistical tests):
		error['the-acronysm'] = '{:.6e}'.format(func.dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all').sum().sum())
	"""

	if args.do_all:
		args.error = ['SDA', 'ADA', 'SSQ', 'CHISQ', 'MNSE', 'PWSD', 'APWSD', 'NPWSD', 'ANPWSD', 'MWUT', 'WMWET', 'TOST', 'DUT']
		"""
		SDA    : Squared Difference of Averages
		ADA    : Absolute Difference of Averages
		SSQ    : Sum of SQuares
		CHISQ  : Chi-Square (Differences divided by data standard deviation)
		MNSE   : Mean Normalized Square Error (Differences divided by data average)
		PWSD   : Pair-Wise Square Deviation
		APWSD  : Absolute Pair-Wise Deviation
		NPWSD  : Normalized Pair-Wise Square Deviation
		ANPWSD : Absolute Normalized Pair-Wise Deviation
		MWUT   : Mann-Whitney U-test (Mann and Whitney, 1947, DOI 10.1214/aoms/1177730491)
		WMWET  : Wellek's Mann-Whitney Equivalence Test (Wellek 1996, DOI 10.1002/bimj.4710380608)
		TOST   : Two one-sided t-tests (Dunnet and Gent, 1977, DOI 10.2307/2529457, as well other authors)
		DUT    : Double Mann-Whitney U-tests (Reviewed in Cornell, 1990, DOI 10.1080/03610929008830433)

		More information in https://pleione.readthedocs.io/en/latest/ObjectiveFunctions.html
		"""

		data_avrg = doavrg(data, len_data)
		data_stdv = dostdv(data, len_data)

		sims_avrg = doavrg(sims, len_sims)
		sims_stdv = dostdv(sims, len_sims)

	# former mean square error, now square difference of means
	if set(args.error).issuperset(set(['SDA'])) or set(args.error).issuperset(set(['MSE'])):
		func = 0

		if not args.do_all:
			data_avrg = doavrg(data, len_data)
			sims_avrg = doavrg(sims, len_sims)

		func = (data_avrg - sims_avrg)**2

		error['SDA'] = '{:.6e}'.format(func.dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all').sum().sum())

	# former mean absolute error, now absolute value of the difference of means
	if set(args.error).issuperset(set(['ADA'])) or set(args.error).issuperset(set(['MAE'])):
		func = 0

		if not args.do_all:
			data_avrg = doavrg(data, len_data)
			sims_avrg = doavrg(sims, len_sims)

		func = abs(data_avrg - sims_avrg)

		error['ADA'] = '{:.6e}'.format(func.dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all').sum().sum())

	# sum of squares (from BioNetFit paper)
	if set(args.error).issuperset(set(['SSQ'])):
		func = 0

		for i in range(len_data):
			for j in range(len_sims):
				func += (data.loc[i] - sims.loc[j])**2

		error['SSQ'] = '{:.6e}'.format(func.dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all').sum().sum())

	# chi-square (from BioNetFit paper)
	if set(args.error).issuperset(set(['CHISQ'])):
		func = 0

		if not args.do_all:
			data_stdv = dostdv(data, len_data)

		for i in range(len_data):
			for j in range(len_sims):
				func += ((data.loc[i] - sims.loc[j]).divide(data_stdv))**2

		error['CHISQ'] = '{:.6e}'.format(func.dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all').sum().sum())

	# mean normalized square error (from BioNetFit paper)
	if set(args.error).issuperset(set(['MNSE'])):
		func = 0

		if not args.do_all:
			data_avrg = doavrg(data, len_data)

		for i in range(len_data):
			for j in range(len_sims):
				func += ((data.loc[i] - sims.loc[j]).divide(data_avrg))**2

		error['MNSE'] = '{:.6e}'.format(func.replace([numpy.inf, -numpy.inf], numpy.nan).dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all').sum().sum())

	# pair-wise square deviation
	if set(args.error).issuperset(set(['PWSD'])):
		func = 0

		for i in range(len_data):
			for j in range(len_sims):
				func += ((data.loc[i] - sims.loc[j])**2).divide(len_data * len_sims)

		error['PWSD'] = '{:.6e}'.format(func.dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all').sum().sum())

	# pair-wise absolute deviation
	if set(args.error).issuperset(set(['APWSD'])):
		func = 0

		for i in range(len_data):
			for j in range(len_sims):
				func += (abs(data.loc[i] - sims.loc[j])).divide(len_data * len_sims)

		error['APWSD'] = '{:.6e}'.format(func.dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all').sum().sum())

	# normalized pair-wise square deviation (also implemented in BioNetFit as equation 3, but not normalized by the number of data * sims)
	if set(args.error).issuperset(set(['NPWSD'])):
		func = 0

		for i in range(len_data):
			for j in range(len_sims):
				func += (((data.loc[i] - sims.loc[j]).divide(data.loc[i]))**2).divide(len_data * len_sims)

		error['NPWSD'] = '{:.6e}'.format(func.replace([numpy.inf, -numpy.inf], numpy.nan).dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all').sum().sum())

	# normalized pair-wise absolute deviation
	if set(args.error).issuperset(set(['ANPWSD'])):
		func = 0

		for i in range(len_data):
			for j in range(len_sims):
				func += (abs((data.loc[i] - sims.loc[j]).divide(data.loc[i]))).divide(len_data * len_sims)

		error['ANPWSD'] = '{:.6e}'.format(func.replace([numpy.inf, -numpy.inf], numpy.nan).dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all').sum().sum())

	"""
	Wellek's Mann-Whitney Equivalence Test.
	Based on mawi.R script from the EQUIVNONINF package
	modifications done to perform the test "vectorized"
	(it compares two matrices; the first has all exp data, the second all the simulations)
	"""
	if set(args.error).issuperset(set(['WMWET'])):
		from scipy.stats import ncx2
		# useful variables (namespace identical to mawi.R script)
		m = len_data # x = data
		n = len_sims # y = sims
		eps1_ = .3129 # Wellek's paper
		eps2_ = .2661 # Wellek's paper
		eqctr = 0.5 + (eps2_ - eps1_)/2
		eqleng = eps1_ + eps2_

		# estimators needed for calculations
		wxy = pandas.DataFrame(index = sims.loc[0].index, columns = sims.loc[0].columns).fillna(0)
		pihxxy = pandas.DataFrame(index = sims.loc[0].index, columns = sims.loc[0].columns).fillna(0)
		pihxyy = pandas.DataFrame(index = sims.loc[0].index, columns = sims.loc[0].columns).fillna(0)
		sigmah = pandas.DataFrame(index = sims.loc[0].index, columns = sims.loc[0].columns).fillna(0)

		# ŷ estimator (wxy in mawi.R)
		# equation 1.2 from Wellek 1996 paper
		# for (i in 1:m) for (j in 1:n) wxy <- wxy + trunc(0.5 * (sign(x[i] - y[j]) + 1))
		for i in range(m):
			for j in range(n):
				diff = (data.loc[i] - sims.loc[j])
				diff = diff.dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all')
				diff = diff.apply(numpy.sign)
				diff = diff + 1
				diff = diff.multiply(0.5)
				diff = diff.apply(numpy.trunc)
				# add to ŷ (wxy in mawi.R)
				wxy += diff

		# yFFG estimator (pihxxy in mawi.R)
		# equation 2.5a from Wellek 1996 paper
		#for (i1 in 1:(m - 1)) for (i2 in (i1 + 1):m) for (j in 1:n) pihxxy <- pihxxy + trunc(0.5 * (sign(min(x[i1], x[i2]) - y[j]) + 1))
		for xi1 in range(m - 1):
			for xi2 in range(xi1 + 1, m):
				for xj in range(n):
					diff = data.loc[xi1].where(data.loc[xi1] < data.loc[xi2], data.loc[xi2]) - sims.loc[xj]
					diff = diff.dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all')
					diff = diff.apply(numpy.sign)
					diff = diff + 1
					diff = diff.multiply(0.5)
					diff = diff.apply(numpy.trunc)
					# add to yFGG (pihxxy in mawi.R)
					pihxxy += diff

		# yFGG estimator (pihxyy in mawi.R)
		# equation 2.5b from Wellek 1996 paper
		# for (i in 1:m) for (j1 in 1:(n - 1)) for (j2 in (j1 + 1):n) pihxyy <- pihxyy + trunc(0.5 * (sign(x[i] - max(y[j1], y[j2])) + 1))
		for xi in range(m):
			for xj1 in range(n - 1):
				for xj2 in range(xj1 + 1, n):
					diff = (data.loc[xi] - sims.loc[xj1].where(sims.loc[xj1] > sims.loc[xj2], sims.loc[xj2]))
					diff = diff.dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all')
					diff = diff.apply(numpy.sign)
					diff = diff + 1
					diff = diff.multiply(0.5)
					diff = diff.apply(numpy.trunc)
					# add to yFGG (pihxyy in mawi.R)
					pihxyy += diff

		# in equation 1.2
		wxy = wxy.divide(m * n)
		# in equation 2.5a, inverse of (m choose 2 = 0.5 * (m-1) * m), then divided by n
		pihxxy = pihxxy.multiply(2).divide(m * (m - 1) * n)
		# in equation 2.5b, inverse of (n choose 2 = 0.5 * (n-1) * n), then divided by m
		pihxyy = pihxyy.multiply(2).divide(n * (n - 1) * m)

		# variance estimator sigmah (same name as in mawi.R)
		# equation 2.6 from Wellek 1996 paper
		# sigmah <- sqrt((wxy - (m + n - 1) * wxy^2 + (m - 1) * pihxxy + (n - 1) * pihxyy)/(m * n))
		sigmah = wxy - (wxy**2).multiply(m + n - 1) + pihxxy.multiply(m - 1) + pihxyy.multiply(n - 1)
		sigmah = sigmah.divide(m * n)
		sigmah = sigmah**0.5

		# critical value
		# right hand of inequality 2.8 from Wellek 1996 paper
		phi = ((eqleng/2)/sigmah)**2
		# crit <- sqrt(qchisq(alpha, 1, (eqleng/2/sigmah)^2))
		# Ca(phi) is the square root of the alpha-th quantile of the chi2-distribution with a single degree of freedom and non-centrality parameter phi square
		crit = pandas.DataFrame(data = ncx2.ppf(0.05, 1, phi), index = sims.loc[0].index, columns = sims.loc[0].columns)**.5

		# compare with Z
		# left hand side of the inequality 2.8 from Wellek 1996 paper
		Z = abs((wxy - eqctr).divide(sigmah))
		z = Z.copy(deep = True)

		"""
		we want to maximize the amount of true alternative hypotheses, so
		we purposely changed the values to use the Wellek's test as an objective function to minimize
		"""
		# test the inequality 2.8 from Wellek 1996 paper
		# the test cannot reject null hypothesis: P[X-Y] < .5 - e1 or P[X-Y] > .5 + e2
		Z[z >= crit] = +1.0
		# the null hypothesis is rejected, therefore .5 - e1 < P[X-Y] < .5 + e2
		Z[z < crit] = +0.0

		if args.report:
			print('wxy estimator:\n', wxy, '\n')
			print('pihxxy estimator:\n', pihxxy, '\n')
			print('pihxyy estimator:\n', pihxyy, '\n')
			print('sigmah estimator:\n', sigmah, '\n')
			print('phi matrix:\n', phi, '\n')
			print('critical values:\n', crit, '\n')
			print('Z estimator: \n', Z, '\n')
			print('Wellek\'s test matrix: a zero means data and simulations are equivalents within the threshold\n', Z)

		error['WMWET'] = '{:.0f}'.format(Z.sum().sum())

	# the same as WMWET, but as identical as the Wellek's paper (look for the heaviside function)
	if set(args.error).issuperset(set(['WMWET_paper'])):
		from scipy.stats import ncx2

		eps1_ = .3129 # Wellek's paper
		eps2_ = .2661 # Wellek's paper
		eqctr = 0.5 + (eps2_ - eps1_)/2
		eqleng = eps1_ + eps2_

		# estimators needed for calculations
		wxy = pandas.DataFrame(index = y.loc[0].index, columns = y.loc[0].columns).fillna(0)
		pihxxy = pandas.DataFrame(index = y.loc[0].index, columns = y.loc[0].columns).fillna(0)
		pihxyy = pandas.DataFrame(index = y.loc[0].index, columns = y.loc[0].columns).fillna(0)
		sigmah = pandas.DataFrame(index = y.loc[0].index, columns = y.loc[0].columns).fillna(0)

		# ŷ estimator (wxy in mawi.R)
		# for (i in 1:m) for (j in 1:n) wxy <- wxy + trunc(0.5 * (sign(x[i] - y[j]) + 1))
		for i in range(m):
			for j in range(n):
				diff = (x.loc[i] - y.loc[j]).dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all')
				wxy += numpy.heaviside(diff, 0)

		# yFFG estimator (pihxxy in mawi.R)
		#for (i1 in 1:(m - 1)) for (i2 in (i1 + 1):m) for (j in 1:n) pihxxy <- pihxxy + trunc(0.5 * (sign(min(x[i1], x[i2]) - y[j]) + 1))
		for xi1 in range(m - 1):
			for xi2 in range(xi1 + 1, m):
				for xj in range(n):
					diff1 = (x.loc[xi1] - y.loc[xj]).dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all')
					diff2 = (x.loc[xi2] - y.loc[xj]).dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all')
					pihxxy += numpy.heaviside(diff1, 0) * numpy.heaviside(diff2, 0)

		# yFGG estimator (pihxyy in mawi.R)
		# for (i in 1:m) for (j1 in 1:(n - 1)) for (j2 in (j1 + 1):n) pihxyy <- pihxyy + trunc(0.5 * (sign(x[i] - max(y[j1], y[j2])) + 1))
		for xi in range(m):
			for xj1 in range(n - 1):
				for xj2 in range(xj1 + 1, n):
					diff1 = (x.loc[xi] - y.loc[xj1]).dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all')
					diff2 = (x.loc[xi] - y.loc[xj2]).dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all')
					pihxyy += numpy.heaviside(diff1, 0) * numpy.heaviside(diff2, 0)

		#
		wxy = wxy.divide(m * n)
		pihxxy = pihxxy.multiply(2).divide(m * (m - 1) * n)
		pihxyy = pihxyy.multiply(2).divide(n * (n - 1) * m)

		# variance estimator sigmah (same name as in mawi.R)
		# sigmah <- sqrt((wxy - (m + n - 1) * wxy^2 + (m - 1) * pihxxy + (n - 1) * pihxyy)/(m * n))
		sigmah = wxy - (wxy**2).multiply(m + n - 1) + pihxxy.multiply(m - 1) + pihxyy.multiply(n - 1)
		sigmah = sigmah.divide(m * n)
		sigmah = sigmah**0.5

		# critical value
		# crit <- sqrt(qchisq(alpha, 1, (eqleng/2/sigmah)^2))
		phi = (eqleng/2/sigmah)**2
		crit = pandas.DataFrame(data = ncx2.ppf(0.05, 1, phi), index = y.loc[0].index, columns = y.loc[0].columns)**.5

		# compare with Z
		Z = abs((wxy - eqctr).divide(sigmah))
		z = Z.copy(deep = True)
		Z[z < crit] = +0.0 # the null hypothesis is rejected, therefore .5 - e1 < P[X-Y] < .5 + e2
		Z[z >= crit] = +1.0 # the test cannot reject the null hypothesis: P[X-Y] < .5 - e1 or P[X-Y] > .5 + e2

		if args.report:
			print('wxy estimator:\n', wxy, '\n')
			print('pihxxy estimator:\n', pihxxy, '\n')
			print('pihxyy estimator:\n', pihxyy, '\n')
			print('sigmah estimator:\n', sigmah, '\n')
			print('phi matrix:\n', phi, '\n')
			print('critical values:\n', crit, '\n')
			print('Z estimator: \n', Z, '\n')
			print('Wellek\'s test matrix: a zero means data and simulations are equivalents within the threshold\n', Z)

		error['WMWET_paper'] = '{:.0f}'.format(Z.sum().sum())

	if set(args.error).issuperset(set(['TOST'])):
		print("WARNING: data and/or simulations not necessarily are normal distributions.")
		print("As a test-bed, we consider data and simulations have unequal standard deviations")
		print("See https://www.statsmodels.org/devel/generated/statsmodels.stats.weightstats.ttost_ind.html for more information")
		from statsmodels.stats.weightstats import ttost_ind

		if not args.do_all:
			data_stdv = dostdv(data, len_data)

		# reshape data and sims to allow calculate the test in a for-loop
		tost_sims = numpy.dstack([sims.loc[x] for x in range(len_sims)])
		# since we operate numpy arrays without labels, we must ensure sims and data indexes and columns have the same order
		index = data.loc[0].index
		columns = data.loc[0].columns
		tost_data = numpy.dstack([data.loc[x].reindex(columns = columns, index = index) for x in range(len_data)])

		p = numpy.zeros((len(data_stdv.index), len(data_stdv.columns)))
		row = 0
		for x, y, lim in zip(tost_sims, tost_data, data_stdv.values):
			for col, _ in enumerate(data_stdv.columns):
				p[row, col] = ttost_ind(x[col], y[col], -lim[col], +lim[col])[0]
			row += 1

		# transform matrix of p-values into a non-rejection DataFrame (if p-value less than 5% -> rejects, but set to zero)
		p = pandas.DataFrame(index = index, columns = columns, data = p)
		P = p.copy(deep = True)
		P[p >= .05] = +1.0
		P[p < .05] = +0.0

		if args.report:
			print('Two one-sided t-tests matrix: a zero means data and simulations are equivalents within one standard deviation threshold\n', P)

		error['TOST'] = '{:.0f}'.format(P.sum().sum())

	# Mann-Whitney U-test
	def mwut(data, sims, alternative):
		ucrit = pandas.read_csv(args.crit, sep = None, engine = 'python', header = 0, index_col = 0)
		udata = pandas.DataFrame(index = sims.loc[0].index, columns = sims.loc[0].columns).fillna(0)
		usims = pandas.DataFrame(index = sims.loc[0].index, columns = sims.loc[0].columns).fillna(0)

		for i in range(len_data):
			for j in range(len_sims):
				Diff = (data.loc[i] - sims.loc[j]).dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all')
				diff = Diff.copy(deep = True)
				# transform data
				# if data < sims, count -1.0
				Diff[diff < 0] = -1.0
				# if data > sims, count +1.0
				Diff[diff > 0] = +1.0
				# if data = sims, count +0.5
				Diff[diff == 0] = +0.5
				# count how many times is data < sims (udata and usims are complementary)
				diff = Diff.copy(deep = True)
				udata += Diff[diff == -1.0].fillna(0).divide(-1) + Diff[diff == +0.5].fillna(0)
				usims += Diff[diff == +1.0].fillna(0).divide(+1) + Diff[diff == +0.5].fillna(0)

		if alternative == 'two-sided':
			# bigU is max(udata, usims), where udata and usims are DataFrames
			bigU = udata.where(udata >= usims).fillna(usims.where(usims >= udata))
		if alternative == 'less':
			bigU = udata
		if alternative == 'greater':
			bigU = usims

		U = len_data * len_sims - bigU
		u = U.copy(deep = True)
		# U is significant if it is less than or equal to a critical value
		U[u <= ucrit.loc[len_sims, str(len_data)]] = +1.0
		U[u > ucrit.loc[len_sims, str(len_data)]] = +0.0

		if args.report:
			print('U-estimator for data\n', udata, '\n')
			print('U-estimator for sims\n', usims, '\n')
			if alternative == 'two-sided':
				print('U-test matrix: A one means data and sims are differents\n', U, '\n')
			if alternative == 'less':
				print('U-test matrix: A one means data is smaller than sims (shifted to the right)\n', U, '\n')
			if alternative == 'greater':
				print('U-test matrix: A one means data is greater than sims (shifted to the left)\n', U, '\n')

		return '{:.0f}'.format(U.sum().sum()), U

	if set(args.error).issuperset(set(['MWUT'])):
		if (len_data >= 3 and len_sims >= 3):
			error['MWUT'] = mwut(data, sims, 'two-sided')[0]
		else:
			error['MWUT'] = str(numpy.nan)

	if set(args.error).issuperset(set(['DUT'])):
		if (len_data >= 3 and len_sims >= 3):
			# set what the user wants
			if args.lower is not None and args.upper is None:
				args.upper = args.lower # symmetric equivalence interval
			if args.lower is None and args.upper is not None:
				args.lower = args.upper # symmetric equivalence interval

			if args.lower is None and args.upper is None:
				if not args.do_all:
					if args.stdv == 'sims':
						lower = upper = dostdv(sims, len_sims)
					else:
						lower = upper = dostdv(data, len_data)
				else:
					if args.stdv == 'sims':
						lower = upper = sims_stdv
					else:
						lower = upper = data_stdv

			# divide by factor
			lower = lower / float(args.factor)
			upper = upper / float(args.factor)

			# copy simulations to a temporary variable
			tmp = sims

			# test lower limit
			new_sims = []
			for i in range(len_sims):
				new_sims.append(tmp.loc[i] - lower)
			sims = pandas.concat(new_sims, keys = range(len_sims))

			# test data > sims - lower with one-tail U-test
			LB = mwut(data, sims, 'greater')[1]

			# test upper limit
			new_sims = []
			for i in range(len_sims):
				new_sims.append(tmp.loc[i] + upper)
			sims = pandas.concat(new_sims, keys = range(len_sims))

			# test data < sims + upper with one
			UB = mwut(data, sims, 'less')[1]

			# rejection DataFrame (U-test report with ones true alternative hypotheses)
			# both one-sided tests should reject the null hypotheses
			U = LB * UB
			# However, we minimize the number of non-rejected null hypotheses
			# transform U into a non-rejection DataFrame.
			U = numpy.logical_xor(U.values, 1).astype(int)
			U = pandas.DataFrame(index = LB.index, columns = LB.columns, data = U)

			if args.report:
				print('Double U-test matrix: 1.0 means data and sims are not equivalents if sims are shifted:\n', U, '\n')

			error['DUT'] = '{:.0f}'.format(U.sum().sum())

		else:
			error['DUT'] = str(numpy.nan)
