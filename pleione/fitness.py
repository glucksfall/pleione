# -*- coding: utf-8 -*-

'''
Project "Genetic Algorithm for Rule-Based Models", Rodrigo Santibáñez, 2017 @ Dlab, FCV (rsantibanez@dlab.cl)
A Python implementation of Alberto Martin's Genetic Algorithm, 2016 @ Dlab, FCV (ajmm@dlab.cl)
To be used with KaSim. Please refer to other subprojects for other stochastic simulators support
Citation:
'''

__author__  = 'Rodrigo Santibáñez'
__license__ = 'gpl-3.0'

#import os, numpy, pandas
import numpy, pandas
from scipy.stats import ncx2

def do(args, sims, len_sims, data, len_data, error):
	"""
	# Fitness Calculation Template:
	if set(args.error).issuperset(set(['the-acronysm'])):
		func = 0
		func = an algebraic expression combining the data average (data_avrg), data variance (data_stdv), simulation average (sims_stdv),
		single experimental files (data.loc[i]) and/or simulation files (sims.loc[i]).
		# Please consider this variables are DataFrames, meaning that division is a method (pandas.DataFrame.division)
		# Please consider use data.loc[i] and sims.loc[i] if average or standard deviation values are needed from them (as in SDM)
		# drop NaN values (from experimental data without simulation point or vice-versa), sum the two dimensions, and return a 6 float points scientific notation number
		error['acronysm'] = '{:.6e}'.format(func.dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all').sum().sum())
	"""

	# former mean square error, now square difference of means
	if set(args.error).issuperset(set(['SDA'])) or set(args.error).issuperset(set(['MSE'])):
		func = 0

		data_avrg = 0
		for i in range(len_data):
			data_avrg += data.loc[i].divide(len_data)

		sims_avrg = 0
		for j in range(len_sims):
			sims_avrg += sims.loc[j].divide(len_sims)

		func = (data_avrg - sims_avrg)**2
		func = func.dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all').sum().sum()

		error['SDA'] = '{:.6e}'.format(func)

	# former mean absolute error, now absolute value of the difference of means
	if set(args.error).issuperset(set(['ADA'])) or set(args.error).issuperset(set(['MAE'])):
		func = 0

		data_avrg = 0
		for i in range(len_data):
			data_avrg += data.loc[i].divide(len_data)

		sims_avrg = 0
		for j in range(len_sims):
			sims_avrg += sims.loc[j].divide(len_sims)

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

		data_avrg = 0
		for i in range(len_data):
			data_avrg += data.loc[i].divide(len_data)

		data_stdv = 0
		if len_data > 1:
			for i in range(len_data):
				data_stdv += ((data.loc[i] - data_avrg)**2).divide(len_data - 1).replace(0., 1.) # replications with stdv = 0 are problematic
		else:
			data_stdv = stdv**2

		for i in range(len_data):
			for j in range(len_sims):
				func += ((data.loc[i] - sims.loc[j]).divide(data_stdv**0.5))**2

		error['CHISQ'] = '{:.6e}'.format(func.dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all').sum().sum())

	# mean normalized square error (from BioNetFit paper)
	if set(args.error).issuperset(set(['MNSE'])):
		func = 0

		data_avrg = 0
		for i in range(len_data):
			data_avrg += data.loc[i].divide(len_data)

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

	# Mann-Whitney U-test
	if set(args.error).issuperset(set(['MWUT'])):
		if ((len_data >= 3 and len_sims >= 3) or (len_data >= 2 and len_sims >= 5)):

			ucrit = pandas.read_table(args.crit[0], sep = None, engine = 'python', header = 0, index_col = 0)
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

			# U is significant if it is less than or equal to the table value
			U = len_data * len_sims - udata.where(udata >= usims).fillna(usims.where(usims >= udata))
			u = U.copy(deep = True)
			U[u <= ucrit.loc[len_sims, str(len_data)]] = +1.0
			U[u > ucrit.loc[len_sims, str(len_data)]] = +0.0

			if args.report:
				print('U-test matrix: 1.0 means distributions are differents\n', U, '\n')

			error['MWUT'] = '{:.0f}'.format(U.sum().sum())

		else:
			error['MWUT'] = str(numpy.nan)

	# Wellek's Mann-Whitney Equivalence Test.
	# Based on mawi.R script from the EQUIVNONINF package
	if set(args.error).issuperset(set(['WMWET'])):
		# set R HOME and import stats R package to use the qchisq function
		# (because the loc arg from scipy.stats.distributions.chi2.ppf gives weird results)
		#from rpy2.robjects.packages import importr
		#os.environ['R_HOME'] = args.r_path
		#os.environ['LD_LIBRARY_PATH'] = args.r_libs
		#stats = importr('stats')

		# useful variables
		m = len_data
		n = len_sims
		e1 = .3129 # Wellek's paper
		e2 = .2661 # Wellek's paper
		e1 = .1382 # mawi.R script
		e2 = .2602 # mawi.R script
		eqctr = .5 + (e2 - e1)/2
		eqleng = e1 + e2

		# estimators needed for calculations
		y = pandas.DataFrame(index = sims.loc[0].index, columns = sims.loc[0].columns).fillna(0)
		yFFG = pandas.DataFrame(index = sims.loc[0].index, columns = sims.loc[0].columns).fillna(0)
		yFGG = pandas.DataFrame(index = sims.loc[0].index, columns = sims.loc[0].columns).fillna(0)
		sigmah = pandas.DataFrame(index = sims.loc[0].index, columns = sims.loc[0].columns).fillna(0)

		# ŷ estimator
		for i in range(m):
			for j in range(n):
				diff = (data.loc[i] - sims.loc[j])
				diff = diff.dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all')
				diff = diff.apply(numpy.sign)
				diff = diff + 1
				diff = diff.multiply(.5)
				# trunc R function
				diff[diff < 0] = diff.apply(numpy.floor)
				diff[diff > 0] = diff.apply(numpy.ceil)
				# add to ŷ (wxy in mawi.R)
				y += diff

		y = y.divide(m*n)
		if args.report:
			print('wxy estimator:\n', y, '\n')

		# yFGG estimator
		for xi in range(m):
			for xj1 in range(n - 1):
				for xj2 in range(xj1 + 1, n):
					diff = (data.loc[xi] - sims.loc[xj1].where(sims.loc[xj1] > sims.loc[xj2], sims.loc[xj2]))
					diff = diff.dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all')
					diff = diff.apply(numpy.sign)
					diff = diff + 1
					diff = diff.multiply(.5)
					# trunc R function
					diff[diff < 0] = diff.apply(numpy.floor)
					diff[diff > 0] = diff.apply(numpy.ceil)
					# add to yFGG (pihxyy in mawi.R)
					yFGG += diff

		yFGG = (yFGG*2).divide(n*(n-1)*m)
		if args.report:
			print('pihxyy estimator:\n', yFGG, '\n')

		# yFFG estimator
		for xi1 in range(m - 1):
			for xi2 in range(xi1 + 1, m):
				for xj in range(n):
					diff = data.loc[xi1].where(data.loc[xi1] < data.loc[xi2], data.loc[xi2]) - sims.loc[xj]
					diff = diff.apply(numpy.sign)
					diff = diff + 1
					diff = diff.multiply(.5)
					# trunc R function
					diff[diff < 0] = diff.apply(numpy.floor)
					diff[diff > 0] = diff.apply(numpy.ceil)
					# add to yFGG (pihxxy in mawi.R)
					yFFG += diff

		yFFG = (yFFG*2).divide(m*(m-1)*n)
		if args.report:
			print('pihxxy estimator:\n', yFFG, '\n')

		# variance estimator sigmah (same name as mawi.R)
		sigmah = (y - (y**2).multiply(m + n - 1) + yFFG.multiply(m - 1) + yFGG.multiply(n - 1)).divide(m*n)
		sigmah = sigmah**.5
		if args.report:
			print('sigmah estimator:\n', sigmah, '\n')

		# critical value
		crit = []
		phi = (eqleng/2/sigmah)**2
		#print(phi)
		phi = phi.values
		a, b = numpy.shape(phi)
		for value in phi.reshape((a*b, 1)):
			if not numpy.isnan(value[0]) or not numpy.isinf(value[0]):
				#rvalue = stats.qchisq(0.05, 1, float(value[0]))[0]
				#print('rvalue', rvalue)
				pvalue = ncx2.ppf(0.05, 1, float(value[0]))
				#print('ncx2', pvalue)
				crit.append(pvalue)
			else:
				crit.append(numpy.nan)
		crit = numpy.asarray(crit).reshape((a, b))
		crit = pandas.DataFrame(data = crit, index = sims.loc[0].index, columns = sims.loc[0].columns)**.5

		if args.report:
			print('critical values:\n', crit, '\n')

		# compare with Z
		Z = abs((y - eqctr).divide(sigmah))
		z = Z.copy(deep = True)
		#print(Z[Z >= crit])
		#print(Z[Z < crit])

		"""
		we purposely changed the values to minimize the function
		we want to minimize the amount of null hypothesis
		"""
		# the test cannot reject null hypothesis: P[X-Y] < .5 - e1 or P[X-Y] > .5 + e2
		Z[z >= crit] = +1.0
		# the test is rejected, therefore .5 - e1 < P[X-Y] < .5 + e2
		Z[z < crit] = +0.0
		#print(numpy.count_nonzero(numpy.isinf(Z)))

		Z = Z.replace([numpy.inf, -numpy.inf], numpy.nan)
		if args.report:
			print('Wellek\'s test matrix: a zero means distributions are equivalents\n', Z)

		error['WMWET'] = '{:.0f}'.format(Z.sum().sum())
