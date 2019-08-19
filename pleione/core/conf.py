# -*- coding: utf-8 -*-

'''
Project "Genetic Algorithm for Rule-Based Models", Rodrigo Santibáñez, 2017 @ Dlab, FCV (rsantibanez@dlab.cl)
A Genetic Algorithm inspired by Alberto Martin's Genetic Algorithm, 2016 @ Dlab, FCV (ajmm@dlab.cl)
Citation:
'''

__author__  = 'Rodrigo Santibáñez'
__license__ = 'gpl-3.0'

import numpy, re
from pleione.core.random import random

def set_regex(simulator):
	if simulator == 'KaSim v4':
		regex = '%\w+: \'(\w+)\' ' \
		'([-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?)\s+(?:\/\/|#)\s+' \
		'(\w+)\[([-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?)\s+' \
		'([-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?)\]\s+' \
		'(\w+)\[([-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?)\s+' \
		'([-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?)\]' \
		'(?:\s+)?([-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+))?\n'

	return regex

def configurate(opts, **kwargs):
	error_msg = ''

	# read the model
	data = []
	with open(opts['model'], 'r') as infile:
		for line in infile:
			data.append(line)

	# find variables to parameterize
	regex = set_regex(opts['simulator'])

	num_pars = 0
	parameters = {}

	for line in range(len(data)):
		matched = re.match(regex, data[line])
		if matched:
			num_pars += 1
			parameters[line] = [
				'par',
				matched.group(1), # parameter name
				matched.group(2), # original value
				matched.group(3), # initial distribution
				matched.group(4), # lower bound or mean
				matched.group(5), # upper bound or standard deviation
				matched.group(6), # mutation distribution
				matched.group(7), # lower bound or probability
				matched.group(8), # upper bound or factor
				matched.group(9), # specific mutation rate (optional)
				]

			# Check validity of parameters
			if matched.group(3) == 'loguniform':
				if float(matched.group(4)) == 0.0:
					error_msg += 'Lower bound for parameter {:s} initial population cannot be zero.'.format(matched.group(1))
			if matched.group(6) == 'loguniform':
				if float(matched.group(4)) == 0.0:
					error_msg += 'Lower bound for parameter {:s} search space cannot be zero.'.format(matched.group(1))
			if matched.group(6) == 'factor':
				if float(matched.group(7)) > 1.0:
					error_msg += 'Mutation probability for parameter {:s} must be a float between zero and one.'.format(matched.group(1))
				if float(matched.group(8)) > 1.0:
					error_msg += 'Mutation foldchange for parameter {:s} must be a float between zero and one.'.format(matched.group(1))

			if matched.group(9) is not None and float(matched.group(9)) > 1.0:
				error_msg += 'Specific mutation probability for parameter {:s} must be a float between zero and one.'.format(matched.group(1))

		else:
			parameters[line] = data[line]

	if num_pars == 0:
		error_msg += 'No variables to parameterize.\n' \
			'Check if selected variables follow the regex (See Manual).'

	# print error
	if error_msg != '':
		print(error_msg)
		raise ValueError(error_msg)

	return parameters

def populate(opts):
	# slice dictionary
	legacy = opts['legacy']
	parameters = opts['parameters']

	# 'parameters' dictionary stores everything in the model, particularly the parameters to fit
	par_keys = list(parameters.keys())

	population = {}
	model_string = 'model_{:0' + str(len(str(opts['num_iter']))) + 'd}' + '_{:0' + str(len(str(opts['pop_size']))) + 'd}'
	for ind in range(opts['pop_size']):
		population['model', ind] = model_string.format(0, ind)
		population['error', ind] = opts['max_error']

		for line in range(len(par_keys)):
			if parameters[line][0] == 'par':
				lower = mean = float(parameters[par_keys[line]][4])
				upper = stdv = float(parameters[par_keys[line]][5])

				if parameters[par_keys[line]][3] == 'uniform':
					population[line, ind] = random.uniform(legacy, lower, upper)
				elif parameters[par_keys[line]][3] == 'loguniform':
					population[line, ind] = numpy.exp(random.uniform(legacy, numpy.log(lower), numpy.log(upper)))
				elif parameters[par_keys[line]][3] == 'lognormal':
					population[line, ind] = random.lognormal(legacy, mean, stdv)
				else:
					raise ValueError('Use uniform/loguniform/lognormal for a valid range to look for parameter values at the first iteration.')

	return population
