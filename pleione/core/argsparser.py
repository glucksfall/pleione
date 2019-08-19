# -*- coding: utf-8 -*-

'''
Project "Genetic Algorithm for Rule-Based Models", Rodrigo Santibáñez, 2017 @ Dlab, FCV (rsantibanez@dlab.cl)
A Genetic Algorithm inspired by Alberto Martin's Genetic Algorithm, 2016 @ Dlab, FCV (ajmm@dlab.cl)
Citation:
'''

__author__  = 'Rodrigo Santibáñez'
__license__ = 'gpl-3.0'

import argparse, sys

def argsparser(**kwargs):
	parser = argparse.ArgumentParser(description = 'Perform a calibration of a RBM employing a Genetic Algorithm.')

	# required arguments
	parser.add_argument('--model'  , metavar = 'path' , type = str  , required = True , nargs = 1  , help = 'RBM with tagged variables to parameterize')
	if not kwargs['simulator'] == 'BNG2':
		parser.add_argument('--final'  , metavar = 'float', type = str  , required = True , nargs = 1  , help = 'limit time to simulate')
		parser.add_argument('--steps'  , metavar = 'float', type = str  , required = True , nargs = 1  , help = 'time steps to simulate')

	# choose one or more fitness functions
	parser.add_argument('--error'  , metavar = 'str'  , type = str  , required = True , nargs = '+', help = 'list of supported fit functions')
	parser.add_argument('--data'   , metavar = 'path' , type = str  , required = True , nargs = '+', help = 'data files to parameterize')

	# useful paths
	parser.add_argument('--python' , metavar = 'path' , type = str  , required = False, default = '~/bin/python3' , help = 'python path, default ~/bin/python3')
	if kwargs['simulator'] == 'BNG2' or kwargs['simulator'] == 'NFsim':
		parser.add_argument('--bng2'   , metavar = 'path' , type = str  , required = False, default = '~/bin/bng2'    , help = 'BioNetGen path, default ~/bin/bng2')
		parser.add_argument('--nfsim'  , metavar = 'path' , type = str  , required = False, default = '~/bin/nfsim'   , help = 'NFsim path, default ~/bin/nfsim')
	if kwargs['simulator'] == 'KaSim v4':
		parser.add_argument('--kasim'  , metavar = 'path' , type = str  , required = False, default = '~/bin/kasim4'  , help = 'KaSim path, default ~/bin/kasim4')
	if kwargs['simulator'] == 'PISKaS':
		parser.add_argument('--piskas' , metavar = 'path' , type = str  , required = False, default = '~/bin/piskas'  , help = 'PISKaS path, default ~/bin/piskas')

	# distribute computation with SLURM, otherwise with python multiprocessing API
	parser.add_argument('--slurm'  , metavar = 'str'  , type = str  , required = False, default = None            , help = 'SLURM partition to use, default None')
	parser.add_argument('--sbatch' , metavar = 'str'  , type = str  , required = False, default = ''              , help = 'explicit configuration for sbatch, e.g. --mem-per-cpu 5G allowing kasim use 5GB RAM per simulation')

	# general options
	if not sys.platform.startswith('linux'):
		parser.add_argument('--seed'   , metavar = 'int'  , type = int  , required = True , default = None            , help = 'random number generator seed, default None')
	else:
		parser.add_argument('--seed'   , metavar = 'int'  , type = int  , required = False, default = None            , help = 'random number generator seed, default None')
	parser.add_argument('--iter'   , metavar = 'int'  , type = int  , required = False, default = 100             , help = 'number of iterations, default 100')
	parser.add_argument('--inds'   , metavar = 'int'  , type = int  , required = False, default = 100             , help = 'number of individuals per iteration, default 100')
	parser.add_argument('--sims'   , metavar = 'int'  , type = int  , required = False, default = 10              , help = 'number of simulations per individual, default 100')
	parser.add_argument('--best'   , metavar = 'int'  , type = int  , required = False, default = 10              , help = 'size of elite individuals, default 10.')
	parser.add_argument('--swap'   , metavar = 'float', type = float, required = False, default = 0.50            , help = 'Q1: global parameter swap probability, default 0.5')
	parser.add_argument('--cross'  , metavar = 'str'  , type = str  , required = False, default = 'multiple'      , help = 'Type of crossover: multiple|single point, default multiple')
	parser.add_argument('--rate'   , metavar = 'float', type = float, required = False, default = 0.50            , help = 'Q2: global parameter mutation probability, default 0.5')
	parser.add_argument('--dist'   , metavar = 'str'  , type = str  , required = False, default = 'inverse'       , help = 'parent selection inverse|uniform, default inverse')
	parser.add_argument('--self'   , metavar = 'False', type = str  , required = False, default = False           , help = 'self recombination True|False, default False')
	parser.add_argument('--crit'   , metavar = 'path' , type = str  , required = False, default = None            , help = 'table of Mann-Whitney U-test critical values, default None')
	parser.add_argument('--prec'   , metavar = 'str'  , type = str  , required = False, default = '7g'            , help = 'precision and format of parameter values, default 7g')

	# other options
	if kwargs['simulator'] == 'KaSim v4':
		parser.add_argument('--syntax' , metavar = 'str'  , type = str  , required = False, default = '4'             , help = 'KaSim syntax, default 4')
	if kwargs['simulator'] == 'KaSim v4-beta':
		parser.add_argument('--binary' , metavar = 'str'  , type = str  , required = False, default = 'model'         , help = 'KaSim binary prefix, default model')
	if kwargs['simulator'] == 'NFsim':
		parser.add_argument('--equil'  , metavar = 'float', type = float, required = False, default = 0               , help = 'equilibrate model before running the simulation, default 0')
	if kwargs['simulator'] == 'PISKaS':
		parser.add_argument('--sync'   , metavar = 'float', type = str  , required = False, default = '1.0'           , help = 'time period to syncronize compartments, default 1.0')
	parser.add_argument('--output' , metavar = 'str'  , type = str  , required = False, default = 'outmodels'     , help = 'ranking files prefixes, default outmodels')
	parser.add_argument('--results', metavar = 'str'  , type = str  , required = False, default = 'results'       , help = 'output folder where to move the results, default results')
	parser.add_argument('--parsets', metavar = 'str'  , type = str  , required = False, default = 'individuals'   , help = 'folder to save the generated models, default individuals')
	parser.add_argument('--rawdata', metavar = 'str'  , type = str  , required = False, default = 'simulations'   , help = 'folder to save the simulations, default simulations')
	parser.add_argument('--fitness', metavar = 'str'  , type = str  , required = False, default = 'goodness'      , help = 'folder to save the goodness of fit, default goodness')
	parser.add_argument('--ranking', metavar = 'str'  , type = str  , required = False, default = 'ranking'       , help = 'folder to save the ranking summaries, default ranking')

	# TO BE DEPRECATED, only with publishing purposes.
	# the random standard library does not have a random.choice with an optional probability list, therefore, Pleione uses numpy.random.choice
	parser.add_argument('--legacy' , metavar = 'True' , type = str  , required = False, default = None            , help = 'use random.random instead of the default numpy.random library')
	# If the user wants to know the behavior of other functions, the option --dev should be maintained
	parser.add_argument('--dev'    , metavar = 'True' , type = str  , required = False, default = None            , help = 'calculate all fitness functions True|False, default False')

	args = parser.parse_args()

	# check some options
	if args.crit is None:
		if set(args.error).issuperset(set(['MWUT'])):
			parser.error('--error MWUT requires --crit file')
		if args.dev:
			parser.error('--dev requires --crit file')
		args.crit = 'dummy-file.txt' # the file is not read by the error calculation script

	if args.seed is None:
		args.seed = int.from_bytes(os.urandom(4), byteorder = 'big')

	if args.legacy is not None and args.dist == 'inverse':
		parser.error('legacy uses the random standard library that don\'t support a non-uniform random choice.\n' \
			'Please delete --legacy True or set it to False.')

	return args
