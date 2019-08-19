# -*- coding: utf-8 -*-

'''
Project "Genetic Algorithm for Rule-Based Models", Rodrigo Santibáñez, 2017 @ Dlab, FCV (rsantibanez@dlab.cl)
A Genetic Algorithm inspired by Alberto Martin's Genetic Algorithm, 2016 @ Dlab, FCV (ajmm@dlab.cl)
Citation:
'''

__author__  = 'Rodrigo Santibáñez'
__license__ = 'gpl-3.0'

# core functions
from pleione.core.argsparser import argsparser
from pleione.core.opts import ga_opts
from pleione.core.random import random
from pleione.core.checks import safe_checks
from pleione.core.utils import clean, backup
from pleione.core.conf import configurate, populate

# genetic algorithm
from pleione.ga.steps import simulate, evaluate, ranking, mutate

def main(**kwargs):
	# general options
	args = argsparser(**kwargs)
	opts = ga_opts(args, **kwargs)
	seed = random.seed(opts['seed'], opts['legacy'])

	# perform checks prior to any calculation
	safe_checks(opts)

	# clean the working directory
	clean(opts)

	# read model configuration
	opts['parameters'] = configurate(opts)
	# generate first population
	opts['population'] = populate(opts)

	# main Genetic Algorithm
	for iter in range(opts['num_iter']):
		opts['iter'] = iter
		opts['population'] = simulate(opts)
		opts['population'] = evaluate(opts)
		opts['population'] = ranking(opts)
		opts['population'] = mutate(opts)

	# move and organize results into subfolders
	backup(opts)
