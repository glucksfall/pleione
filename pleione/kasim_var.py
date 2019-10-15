# -*- coding: utf-8 -*-

'''
Project "Genetic Algorithm for Rule-Based Models", Rodrigo Santibáñez, 2017 @ Dlab, FCV (rsantibanez@dlab.cl)
A Python implementation of Alberto Martin's Genetic Algorithm, 2016 @ Dlab, FCV (ajmm@dlab.cl)
To be used with KaSim. Please refer to other subprojects for other stochastic simulators support
Citation: Pleione: A tool for statistical and multi-objective calibration of Rule-based models. Scientific Reports (2019)
DOI:
'''

__author__  = 'Rodrigo Santibáñez'
__license__ = 'gpl-3.0'
__software__ = 'kasim-v4.0'

import argparse, glob, multiprocessing, os, random, re, shutil, subprocess, sys, time
import pandas, numpy

class custom:
	class random:
		def seed(number):
			if args.legacy:
				random.seed(args.seed)
			else:
				numpy.random.seed(args.seed)

		def random():
			if args.legacy:
				return random.random()
			else:
				return numpy.random.random()

		def uniform(lower, upper):
			if args.legacy:
				return random.uniform(lower, upper)
			else:
				return numpy.random.uniform(lower, upper, None)

		def lognormal(lower, upper):
			if args.legacy:
				return random.lognormvariate(lower, upper)
			else:
				return numpy.random.lognormal(lower, upper, None)

def safe_checks():
	error_msg = ''
	if shutil.which(opts['python']) is None:
		error_msg += 'python3 (at {:s}) can\'t be called to perform error calculation.\n' \
			'You could use --python {:s}\n'.format(opts['python'], shutil.which('python3'))

	# check for simulators
	#if shutil.which(opts['bng2']) is None:
		#error_msg += 'BNG2 (at {:s}) can\'t be called to perform simulations.\n' \
			#'Check the path to BNG2.'.format(opts['bng2'])
	if shutil.which(opts['kasim']) is None:
		error_msg += 'KaSim (at {:s}) can\'t be called to perform simulations.\n' \
			'Check the path to KaSim.'.format(opts['kasim'])
	#if shutil.which(opts['nfsim']) is None:
		#error_msg += 'NFsim (at {:s}) can\'t be called to perform simulations.\n' \
			#'Check the path to NFsim.'.format(opts['nfsim'])
	#if shutil.which(opts['piskas']) is None:
		#error_msg += 'PISKaS (at {:s}) can\'t be called to perform simulations.\n' \
			#'Check the path to PISKaS.'.format(opts['piskas'])

	# check for slurm
	if opts['slurm'] is not None or opts['slurm'] == '':
		if not sys.platform.startswith('linux'):
			error_msg += 'SLURM do not support WindowsOS and macOS (https://slurm.schedmd.com/platforms.html)\n'
		else:
			if shutil.which('sinfo') is None:
				error_msg += 'You specified a SLURM partition but SLURM isn\'t installed on your system.\n' \
					'Delete --slurm to use the python multiprocessing API or install SLURM (https://pleione.readthedocs.io/en/latest/SLURM.html)\n'
			else:
				cmd = 'sinfo -hp {:s}'.format(opts['slurm'])
				cmd = re.findall(r'(?:[^\s,"]|"+(?:=|\\.|[^"])*"+)+', cmd)
				out, err = subprocess.Popen(cmd, shell = False, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
				if out == b'':
					error_msg += 'You specified an invalid SLURM partition.\n' \
						'Please, use --slurm $SLURM_JOB_PARTITION or delete --slurm to use the python multiprocessing API.\n'

	# check if model file exists
	if not os.path.isfile(opts['model']):
		error_msg += 'The "{:s}" file cannot be opened.\n' \
			'Please, check the path to the model file.\n'.format(opts['model'])

	# check if data files exist
	if len(opts['data']) == 1: # shlex
		if len(glob.glob(opts['data'][0])) == 0:
			error_msg += 'The path "{:s}" is empty.\n' \
				'Please, check the path to the data files.\n'.format(opts['data'][0])
	else:
		for data in opts['data']: # each file was declared explicitly
			if not os.path.isfile(data):
				error_msg += 'The "{:s}" file cannot be opened.\n' \
					'Please, check the path to the data file.\n'.format(data)

	# check GA options
	if opts['mut_swap'] > 1.0:
		error_msg += 'Parameter swap (or recombination) probability must be a float between zero and one.\n'

	if opts['mut_rate'] > 1.0:
		error_msg += 'Parameter mutation probability must be a float between zero and one.\n'

	# print error
	if error_msg != '':
		print(error_msg)
		raise ValueError(error_msg)

	return 0

def parallelize(cmd):
	proc = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
	out, err = proc.communicate()
	proc.wait()
	return 0

def argsparser():
	parser = argparse.ArgumentParser(description = 'Perform a calibration of a RBM employing a Genetic Algorithm.', \
		epilog = 'cite "Pleione: A tool for statistical and multi-objective calibration of Rule-based models. Scientific Reports (2019)"',
		formatter_class = argparse.RawTextHelpFormatter)

	# required arguments
	parser.add_argument('--model'  , metavar = 'str'  , type = str  , required = True , nargs = 1  , help = 'RBM with tagged variables to parameterize')
	parser.add_argument('--final'  , metavar = 'float', type = str  , required = True , nargs = 1  , help = 'limit time to simulate')
	parser.add_argument('--steps'  , metavar = 'float', type = str  , required = True , nargs = 1  , help = 'time steps to simulate')
	# choose one or more fitness functions
	parser.add_argument('--error'  , metavar = 'str'  , type = str  , required = True , nargs = '+', help = 'list of supported fit functions')
	parser.add_argument('--data'   , metavar = 'str'  , type = str  , required = True , nargs = '+', help = 'data files to parameterize')

	# useful paths
	#parser.add_argument('--bng2'   , metavar = 'path' , type = str  , required = False, default = '~/bin/bng2'    , help = 'BioNetGen path, default ~/bin/bng2')
	parser.add_argument('--kasim'  , metavar = 'path' , type = str  , required = False, default = '~/bin/kasim4'  , help = 'KaSim path, default ~/bin/kasim4')
	#parser.add_argument('--nfsim'  , metavar = 'path' , type = str  , required = False, default = '~/bin/nfsim'   , help = 'NFsim path, default ~/bin/nfsim')
	#parser.add_argument('--piskas' , metavar = 'path' , type = str  , required = False, default = '~/bin/piskas'  , help = 'PISKaS path, default ~/bin/piskas')
	parser.add_argument('--python' , metavar = 'path' , type = str  , required = False, default = '~/bin/python3' , help = 'python path, default ~/bin/python3')

	# distribute computation with SLURM, otherwise with python multiprocessing API
	parser.add_argument('--slurm'  , metavar = 'str'  , type = str  , required = False, default = None            , help = 'SLURM partition to use, default None')

	# general options
	parser.add_argument('--seed'   , metavar = 'int'  , type = int  , required = False, default = None            , help = 'random number generator seed, default None')
	parser.add_argument('--iter'   , metavar = 'int'  , type = int  , required = False, default = 100             , help = 'number of iterations, default 100')
	parser.add_argument('--inds'   , metavar = 'int'  , type = int  , required = False, default = 100             , help = 'number of individuals per iteration, default 100')
	parser.add_argument('--sims'   , metavar = 'int'  , type = int  , required = False, default = 10              , help = 'number of simulations per individual, default 100')
	parser.add_argument('--best'   , metavar = 'int'  , type = int  , required = False, default = 10              , help = 'size of elite individuals, default 10.')
	parser.add_argument('--swap'   , metavar = 'float', type = float, required = False, default = 0.50            , help = 'Q1: global parameter swap probability, default 0.5')
	parser.add_argument('--cross'  , metavar = 'str'  , type = str  , required = False, default = 'multiple'      , help = 'Type of crossover: multiple or single point, default multiple')
	parser.add_argument('--rate'   , metavar = 'float', type = float, required = False, default = 0.50            , help = 'Q2: global parameter mutation probability, default 0.5')
	parser.add_argument('--dist'   , metavar = 'str'  , type = str  , required = False, default = 'inverse'       , help = 'parent selection inverse|uniform, default inverse')
	parser.add_argument('--self'   , metavar = 'False', type = str  , required = False, default = False           , help = 'self recombination True|False, default False')
	parser.add_argument('--crit'   , metavar = 'path' , type = str  , required = False, default = None            , help = 'table of Mann-Whitney U-test critical values, default None')
	parser.add_argument('--prec'   , metavar = 'str'  , type = str  , required = False, default = '7g'            , help = 'precision and format of parameter values, default 7g')

	# other options
	parser.add_argument('--syntax' , metavar = 'str'  , type = str  , required = False, default = '4'             , help = 'KaSim syntax, default 4')
	parser.add_argument('--binary' , metavar = 'str'  , type = str  , required = False, default = 'model'         , help = 'KaSim binary prefix, default model')
	#parser.add_argument('--equil'  , metavar = 'float', type = float, required = False, default = 0               , help = 'equilibrate model before running the simulation, default 0')
	#parser.add_argument('--sync'   , metavar = 'float', type = str  , required = False, default = '1.0'           , help = 'time period to syncronize compartments, default 1.0')
	parser.add_argument('--output' , metavar = 'str'  , type = str  , required = False, default = 'outmodels'     , help = 'ranking files prefixes, default outmodels')
	parser.add_argument('--results', metavar = 'str'  , type = str  , required = False, default = 'results'       , help = 'output folder where to move the results, default results')
	parser.add_argument('--parsets', metavar = 'str'  , type = str  , required = False, default = 'individuals'   , help = 'folder to save the generated models, default individuals')
	parser.add_argument('--rawdata', metavar = 'str'  , type = str  , required = False, default = 'simulations'   , help = 'folder to save the simulations, default simulations')
	parser.add_argument('--fitness', metavar = 'str'  , type = str  , required = False, default = 'goodness'      , help = 'folder to save the goodness of fit, default goodness')
	parser.add_argument('--ranking', metavar = 'str'  , type = str  , required = False, default = 'ranking'       , help = 'folder to save the ranking summaries, default ranking')

	# TO BE DEPRECATED, only with publishing purposes.
	# the random standard library does not have a random.choice with an optional probability list, therefore, Pleione uses numpy.random.choice
	parser.add_argument('--legacy' , metavar = 'True' , type = str  , required = False, default = False           , help = 'use True: random.random instead of False: numpy.random')
	# If the user wants to know the behavior of other functions, the option --dev should be maintained
	parser.add_argument('--dev'    , metavar = 'True' , type = str  , required = False, default = False           , help = 'calculate all fitness functions, default False')

	args = parser.parse_args()

	if args.crit is None:
		if set(args.error).issuperset(set(['MWUT'])):
			parser.error('--error MWUT requires --crit file')
		if args.dev:
			parser.error('--dev requires --crit file')
		args.crit = 'dummy-file.txt' # the file is not read by the error calculation script

	if args.seed is None:
		if sys.platform.startswith('linux'):
			args.seed = int.from_bytes(os.urandom(4), byteorder = 'big')
		else:
			parser.error('pleione requires --seed integer')

	if not args.legacy and args.dist == 'inverse':
		parser.error('legacy uses the random standard library that don\'t support a non-uniform random choice.\n' \
			'Please delete legacy or set to False.')

	return args

def ga_opts():
	return {
		# user defined options
		'model'     : args.model[0],
		'final'     : args.final[0], # not bng2
		'steps'     : args.steps[0], # not bng2
		'error'     : args.error,
		'data'      : args.data,
		#'bng2'      : os.path.expanduser(args.bng2), # bng2, nfsim only
		'kasim'     : os.path.expanduser(args.kasim), # kasim only
		#'piskas'    : os.path.expanduser(args.piskas), # piskas only
		#'nfsim'     : os.path.expanduser(args.nfsim), # nfsim only
		'python'    : os.path.expanduser(args.python),
		'slurm'     : args.slurm,
		'rng_seed'  : args.seed,
		'num_iter'  : args.iter,
		'pop_size'  : args.inds,
		'num_sims'  : args.sims,
		'pop_best'  : args.best,
		'mut_swap'  : args.swap,
		'mut_rate'  : args.rate,
		'dist_type' : args.dist,
		'self_rec'  : args.self,
		'xpoints'   : args.cross,
		'crit_vals' : args.crit,
		'par_fmt'   : args.prec,
		'syntax'    : args.syntax, # kasim only
		'binary'    : args.binary, # kasim only
		#'equil'     : args.equil, # nfsim only
		#'sync'      : args.sync, # piskas only
		'outfile'   : args.output,
		'results'   : args.results,
		'parsets'   : args.parsets,
		'rawdata'   : args.rawdata,
		'fitness'   : args.fitness,
		'ranking'   : args.ranking,
		# non-user defined options
		'home'      : os.getcwd(),
		'null'      : '/dev/null',
		'max_error' : numpy.nan,
		'bin_file'  : args.model[0].split('.')[0] + '.bin', # kasim4 only
		'systime'   : str(time.time()).split('.')[0],
		}

def configurate():
	# read the model
	data = []
	with open(opts['model'], 'r') as infile:
		for line in infile:
			data.append(line)

	# find variables to parameterize
	regex = '%\w+: \'(\w+)\' ' \
		'([-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?)\s+(?:\/\/|#)\s+' \
		'(\w+)\[([-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?)\s+' \
		'([-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?)\]\s+' \
		'(\w+)\[([-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?)\s+' \
		'([-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+)?)\]' \
		'(?:\s+)?([-+]?(?:(?:\d*\.\d+)|(?:\d+\.?))(?:[Ee][+-]?\d+))?\n'

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
					print(error_msg)
					raise ValueError(error_msg)
			if matched.group(6) == 'loguniform':
				if float(matched.group(4)) == 0.0:
					error_msg += 'Lower bound for parameter {:s} search space cannot be zero.'.format(matched.group(1))
					print(error_msg)
					raise ValueError(error_msg)

			if matched.group(6) == 'factor':
				if float(matched.group(7)) > 1.0:
					error_msg += 'Mutation probability for parameter {:s} must be a float between zero and one.'.format(matched.group(1))
					print(error_msg)
					raise ValueError(error_msg)
				if float(matched.group(8)) > 1.0:
					error_msg += 'Mutation foldchange for parameter {:s} must be a float between zero and one.'.format(matched.group(1))
					print(error_msg)
					raise ValueError(error_msg)

			if matched.group(9) is not None and float(matched.group(9)) > 1.0:
				error_msg += 'Specific mutation probability for parameter {:s} must be a float between zero and one.'.format(matched.group(1))
				print(error_msg)
				raise ValueError(error_msg)

		else:
			parameters[line] = data[line]

	if num_pars == 0:
		error_msg += 'No variables to parameterize.\n' \
			'Check if selected variables follow the regex (See Manual).'
		print(error_msg)
		raise ValueError(error_msg)

	return parameters

def populate():
	# 'parameters' dictionary stores everything in the model, particularly the parameters to fit
	par_keys = list(parameters.keys())

	population = {}
	for ind in range(opts['pop_size']):
		population['model', ind] = 'model_000_{:03d}'.format(ind)
		population['error', ind] = opts['max_error']

		for line in range(len(par_keys)):
			if parameters[line][0] == 'par':
				lower = mean = float(parameters[par_keys[line]][4])
				upper = stdv = float(parameters[par_keys[line]][5])

				if parameters[par_keys[line]][3] == 'uniform':
					population[line, ind] = custom.random.uniform(lower, upper)
				elif parameters[par_keys[line]][3] == 'loguniform':
					population[line, ind] = numpy.exp(custom.random.uniform(numpy.log(lower), numpy.log(upper)))
				elif parameters[par_keys[line]][3] == 'lognormal':
					population[line, ind] = custom.random.lognormal(mean, stdv)
				else:
					raise ValueError('Declare uniform/loguniform/lognormal for a valid range to look for parameter values at the first iteration')

	return population

def make_bin():
	job_desc = {
		'nodes'     : 1,
		'ntasks'    : 1,
		'ncpus'     : 1,
		'null'      : opts['null'],
		'partition' : opts['slurm'],
		'job_name'  : 'child_{:s}'.format(opts['systime']),
		'stdout'    : 'stdout_{:s}.txt'.format(opts['systime']),
		'stderr'    : 'stderr_{:s}.txt'.format(opts['systime']),
		'make_bin'  : '{:s} -i {:s} -make-sim {:s} -l 0 -o {:s} -syntax {:s} -mode batch --no-log'.format(opts['kasim'], opts['model'], opts['bin_file'], opts['null'], opts['syntax']),
		}

	# submit to the queue
	squeue = []
	# use SLURM Workload Manager
	if opts['slurm'] is not None:
		cmd = os.path.expanduser('sbatch --no-requeue -p {partition} -N {nodes} -c {ncpus} -n {ntasks} -o {null} -e {null} -J {job_name} --wrap ""{make_bin}""'.format(**job_desc))
		cmd = re.findall(r'(?:[^\s,"]|"+(?:=|\\.|[^"])*"+)+', cmd)
		out, err = subprocess.Popen(cmd, shell = False, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
		while err == sbatch_error:
			out, err = subprocess.Popen(cmd, shell = False, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
		squeue.append(out.decode('utf-8')[20:-1])

	# use multiprocessing.Pool
	else:
		cmd = os.path.expanduser(job_desc['make_bin'])
		cmd = re.findall(r'(?:[^\s,"]|"+(?:=|\\.|[^"])*"+)+', cmd)
		out, err = subprocess.Popen(cmd, shell = False, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()

	# check if job has finished
	while not os.path.exists('{:s}'.format(opts['bin_file'])):
		time.sleep(1)

	return population

def simulate():
	job_desc = {
		'nodes'     : 1,
		'ntasks'    : 1,
		'ncpus'     : 1,
		'null'      : opts['null'],
		'partition' : opts['slurm'],
		'job_name'  : 'child_{:s}'.format(opts['systime']),
		'stdout'    : 'stdout_{:s}.txt'.format(opts['systime']),
		'stderr'    : 'stderr_{:s}.txt'.format(opts['systime']),
		}

	# generate a sh file per model
	par_keys = list(parameters.keys())
	par_string = '-var {:s} {:.' + opts['par_fmt'] + '}\n'
	for ind in range(opts['pop_size']):
		model = population['model', ind]

		if not os.path.exists(model + '.sh'):
			with open(model + '.sh', 'w') as file:
				free_params = ''
				for line in range(len(par_keys)):
					if parameters[line][0] == 'par':
						file.write(par_string.format(parameters[line][1], population[line, ind]))
						free_params += par_string.format(parameters[line][1], population[line, ind]).replace('\n', ' ')

	# submit simulations to the queue
	squeue = []

	for ind in range(opts['pop_size']):
		for sim in range(opts['num_sims']):
			model = population['model', ind]
			output = '{:s}.{:03d}.out.txt'.format(model, sim)

			if not os.path.exists(output):
				# read the parameters (Yeap, I could read the population dictionary, but this was the original implementation)
				#cmd = 'cat {:s}.sh'.format(model)
				#cmd = re.findall(r'(?:[^\s,"]|"+(?:=|\\.|[^"])*"+)+', cmd)
				#out, err = subprocess.Popen(cmd, shell = False, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
				#free_params = out.decode('utf-8').replace('\n', ' ')

				#job_desc['exec_kasim'] = '{:s} -load-sim {:s} -l {:s} -p {:s} -o {:s} -syntax {:s} {:s} -mode batch --no-log'.format( \
					#opts['kasim'], opts['bin_file'], opts['final'], opts['steps'], output, opts['syntax'], free_params)
				job_desc['exec_kasim'] = '{:s} -i {:s} -l {:s} -p {:s} -o {:s} -syntax {:s} {:s} -mode batch --no-log'.format( \
					opts['kasim'], opts['model'], opts['final'], opts['steps'], output, opts['syntax'], free_params)

				# use SLURM Workload Manager
				if opts['slurm'] is not None:
					cmd = os.path.expanduser('sbatch --no-requeue -p {partition} -N {nodes} -c {ncpus} -n {ntasks} -o {null} -e {null} -J {job_name} --wrap ""{exec_kasim}""'.format(**job_desc))
					cmd = re.findall(r'(?:[^\s,"]|"+(?:=|\\.|[^"])*"+)+', cmd)
					out, err = subprocess.Popen(cmd, shell = False, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
					while err == sbatch_error:
						out, err = subprocess.Popen(cmd, shell = False, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
					squeue.append(out.decode('utf-8')[20:-1])

				# use multiprocessing.Pool
				else:
					cmd = os.path.expanduser(job_desc['exec_kasim'])
					cmd = re.findall(r'(?:[^\s,"]|"+(?:=|\\.|[^"])*"+)+', cmd)
					squeue.append(cmd)

	# check if squeued jobs have finished
	if opts['slurm'] is not None:
		for job_id in range(len(squeue)):
			cmd = 'squeue --noheader -j{:s}'.format(squeue[job_id])
			cmd = re.findall(r'(?:[^\s,"]|"+(?:=|\\.|[^"])*"+)+', cmd)
			out, err = subprocess.Popen(cmd, shell = False, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
			while out.count(b'child') > 0 or err == squeue_error:
				time.sleep(1)
				out, err = subprocess.Popen(cmd, shell = False, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()

	# simulate with multiprocessing.Pool
	else:
		with multiprocessing.Pool(multiprocessing.cpu_count() - 1) as pool:
			pool.map(parallelize, sorted(squeue), chunksize = 1)

	return population

def evaluate():
	job_desc = {
		'nodes'     : 1,
		'ntasks'    : 1,
		'ncpus'     : 1,
		'null'      : opts['null'],
		'partition' : opts['slurm'],
		'job_name'  : 'child_{:s}'.format(opts['systime']),
		'stdout'    : 'stdout_{:s}.txt'.format(opts['systime']),
		'stderr'    : 'stderr_{:s}.txt'.format(opts['systime']),
		'doerror'   : '{:s} -m pleione.kasim-doerror --crit {:s}'.format(opts['python'], opts['crit_vals']),
		#'deverror'  : '{:s} -m pleione.kasim-allerror --crit {:s}'.format(opts['python'], opts['crit_vals']),
		}

	# submit error calculations to the queue
	squeue = []

	for ind in range(opts['pop_size']):
		model = population['model', ind]

		data = ' '.join(glob.glob(' '.join(opts['data'])))
		error = ' '.join(opts['error'])
		sims = ' '.join(glob.glob('{:s}.*.out.txt'.format(model)))
		output = '{:s}.txt'.format(model)

		job_desc['calc'] = job_desc['doerror'] + ' --data {:s} --sims {:s} --file {:s} --error {:s}'.format(data, sims, output, error)
		if args.dev:
			job_desc['calc'] = job_desc['doerror'] + ' --data {:s} --sims {:s} --file {:s} --do_all True'.format(data, sims, output)

		# use SLURM Workload Manager
		if opts['slurm'] is not None:
			cmd = 'sbatch --no-requeue -p {partition} -N {nodes} -c {ncpus} -n {ntasks} -o {null} -e {null} -J {job_name} --wrap ""{calc}""'.format(**job_desc)
			cmd = re.findall(r'(?:[^\s,"]|"+(?:=|\\.|[^"])*"+)+', cmd)
			out, err = subprocess.Popen(cmd, shell = False, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
			while err == sbatch_error:
				out, err = subprocess.Popen(cmd, shell = False, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
			squeue.append(out.decode('utf-8')[20:-1])

		# use multiprocessing.Pool
		else:
			cmd = os.path.expanduser(job_desc['calc'])
			cmd = re.findall(r'(?:[^\s,"]|"+(?:=|\\.|[^"])*"+)+', cmd)
			squeue.append(cmd)

	# check if squeued jobs have finished
	if opts['slurm'] is not None:
		for job_id in range(len(squeue)):
			cmd = 'squeue --noheader -j{:s}'.format(squeue[job_id])
			cmd = re.findall(r'(?:[^\s,"]|"+(?:=|\\.|[^"])*"+)+', cmd)
			out, err = subprocess.Popen(cmd, shell = False, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
			while out.count(b'child') > 0 or err == squeue_error:
				time.sleep(1)
				out, err = subprocess.Popen(cmd, shell = False, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()

	# calc error with multiprocessing.Pool
	else:
		with multiprocessing.Pool(multiprocessing.cpu_count() - 1) as pool:
			pool.map(parallelize, sorted(squeue), chunksize = 1)

	return population

def ranking():
	for ind in range(opts['pop_size']):
		with open('{:s}.txt'.format(population['model', ind]), 'r') as file:
			data = pandas.read_csv(file, delimiter = '\t', header = None).set_index(0, drop = False).rename_axis(None, axis = 0).drop(0, axis = 1).rename(columns = {1: 'value'})

		if args.dev:
			fitfunc = list(data.index)
		else:
			fitfunc = opts['error']

		# and store the error in the population dictionary
		for name in range(len(fitfunc)):
			population['{:s}'.format(fitfunc[name]), ind] = data.loc['{:s}'.format(fitfunc[name]), 'value']

	# now that everything is stored in the population dict, we proceed to rank models by the selected error function(s)
	jobs = {}
	rank = {}
	fitfunc = opts['error']

	for name in range(len(fitfunc)):
		for ind in range(opts['pop_size']):
			jobs[population['model', ind]] = population['{:s}'.format(fitfunc[name]), ind]
		rank['{:s}'.format(fitfunc[name])] = sorted(jobs, key = jobs.get, reverse = False)

	for name in range(len(fitfunc)):
		for ind in range(opts['pop_size']):
			jobs[population['model', ind]] += { key : value for value, key in enumerate(rank['{:s}'.format(fitfunc[name])]) }[population['model', ind]]

	# create an 'ordered' list of individuals from the 'population' dictionary by increasing fitness
	rank = sorted(jobs, key = jobs.get, reverse = False)

	# find the index that match best individual with the 'population' dictionary keys and store the rank (a list) in the population dictionary
	ranked_population = []
	for best in range(opts['pop_size']):
		for ind in range(opts['pop_size']):
			if population['model', ind] == rank[best]:
				ranked_population.append(ind)
				break

	population['rank'] = ranked_population

	# save the population dictionary as a report file
	par_keys = list(set([x[0] for x in population.keys() if str(x[0]).isdigit()]))

	if args.dev:
		fitfunc = sorted(list(data.index))

	with open('{:s}_{:03d}.txt'.format(opts['outfile'], iter), 'w') as file:
		file.write('# Output of {:s} {:s}\n'.format(opts['python'], subprocess.list2cmdline(sys.argv[0:])))
		file.write('Elapsed time: {:.0f} seconds\n'.format(time.time() - float(opts['systime'])))
		file.write('iteration: {:03d}\t'.format(iter))
		file.write('error: {:s}\t'.format(','.join(opts['error'])))
		file.write('seed: {:d}\n\n'.format(opts['rng_seed']))

		# header
		file.write('MODEL_ID\t')
		for i in range(len(fitfunc)):
			file.write('{:s}\t'.format(fitfunc[i]))
		for key in range(len(par_keys)):
			file.write('{:s}\t'.format(parameters[par_keys[key]][1].strip()))
		file.write('\n')

		for ind in ranked_population:
			file.write('{:s}\t'.format(population['model', ind]))
			for i in range(len(fitfunc)):
				if fitfunc[i] not in ['MWUT', 'WMWET', 'TOST']:
					file.write('{:.6e}\t'.format(float(population['{:s}'.format(fitfunc[i]), ind])))
				else:
					file.write('{:.0f}\t'.format(float(population['{:s}'.format(fitfunc[i]), ind])))
			for key in range(len(par_keys)):
				file.write('{:.6e}\t'.format(float(population[par_keys[key], ind])))
			file.write('\n')

	return population

def mutate():
	# par_keys stores parameter values only
	par_keys = list(set([key[0] for key in population.keys() if str(key[0]).isdigit()]))

	# slice the population dictionary retrieving the best models, if needed
	ranked_population = population['rank']

	if opts['pop_best'] == 0:
		best_population = population
	else:
		best_population = {}
		for best in range(opts['pop_best']):
			for key in range(len(par_keys)):
				best_population[par_keys[key], best] = population[par_keys[key], ranked_population[best]]
			best_population['model', best] = population['model', ranked_population[best]]
			best_population['error', best] = population['error', ranked_population[best]]

	# fill the population dictionary with the elite, because population is a global variable
	for ind in range(opts['pop_best']):
		for par in range(len(par_keys)):
			population[par_keys[par], ind] = best_population[par_keys[par], ind]
		population['model', ind] = best_population['model', ind]
		population['error', ind] = best_population['error', ind]

	# User defined best population
	top = opts['pop_best']
	if opts['pop_best'] == 0:
		top = opts['pop_size']
	elif opts['pop_best'] == 1:
		opts['self_rec'] == True # allow self recombination to generate descendants from only one parent

	# probability distribution to select parents
	if opts['dist_type'] == 'uniform':
		# Define a uniform probability distribution according to the best population size
		dist = [1 for n in range(1, top + 1)]
	elif opts['dist_type'] == 'inverse':
		# Define an inverse probability distribution according to the best population size
		dist = [1/n for n in range(1, top + 1)]
	prob = numpy.divide(dist, float(numpy.sum(dist)))

	# fill the population dictionary with individuals from the best parents
	for ind in range(opts['pop_best'], opts['pop_size'], 2):
		if opts['pop_best'] == 0:
			# choose two random individuals from the ranked population (index start at zero)
			n1 = numpy.random.choice(ranked_population[0:top], p = prob[0:top])
			n2 = numpy.random.choice(ranked_population[0:top], p = prob[0:top])
			if opts['self_rec'] == False and not opts['pop_size'] == 1:
				while n2 == n1:
					n2 = numpy.random.choice(ranked_population[0:top], p = prob[0:top])

		elif opts['pop_best'] != 0:
			if not args.legacy:
				# choose two random individuals from the best population (reindexed from 0 to 'pop_best' size)
				n1 = numpy.random.choice(range(top), p = prob[0:top])
				n2 = numpy.random.choice(range(top), p = prob[0:top])
			else:
				n1 = random.choice(range(top))
				n2 = random.choice(range(top))

			if opts['self_rec'] == False and not opts['pop_size'] == 1:
				while n2 == n1:
					if args.legacy:
						n2 = numpy.random.choice(range(top), p = prob[0:top])
					else:
						n2 = random.choice(range(top))

		# perform multiple or single crossover
		if opts['xpoints'] == 'multiple':
			for par in range(len(par_keys)):
				# create children
				population[par_keys[par], ind] = best_population[par_keys[par], n1]
				population[par_keys[par], ind + 1] = best_population[par_keys[par], n2]

				# swap parameter values using a probability threshold
				if opts['mut_swap'] >= custom.random.random():
					population[par_keys[par], ind] = best_population[par_keys[par], n2]
					population[par_keys[par], ind + 1] = best_population[par_keys[par], n1]

		elif opts['xpoints'] == 'single':
			point = custom.random.uniform(0, len(par_keys))
			for par in range(len(par_keys)):
				# create children and do not swap parameters!
				if par <= point:
					population[par_keys[par], ind] = best_population[par_keys[par], n1]
					population[par_keys[par], ind + 1] = best_population[par_keys[par], n2]
				else:
					population[par_keys[par], ind] = best_population[par_keys[par], n2]
					population[par_keys[par], ind + 1] = best_population[par_keys[par], n1]

		# include the model id
		population['model', ind] = 'model_{:03d}_{:03d}'.format(iter + 1, ind)
		population['model', ind + 1] = 'model_{:03d}_{:03d}'.format(iter + 1, ind + 1)

		# include the error in the population dictionary
		population['error', ind] = opts['max_error']
		population['error', ind + 1] = opts['max_error']

	# mutate parameter values
	for ind in range(opts['pop_best'], opts['pop_size']):
		for par in range(len(par_keys)):
			if parameters[par_keys[par]][6] == 'factor':
				if float(parameters[par_keys[par]][7]) >= custom.random.random():
					population[par_keys[par], ind] = population[par_keys[par], ind] * \
						custom.random.uniform(1.0 - float(parameters[par_keys[par]][8]), 1.0 + float(parameters[par_keys[par]][8]))

			elif parameters[par_keys[par]][6] == 'uniform' or parameters[par_keys[par]][6] == 'loguniform':
				lower = float(parameters[par_keys[par]][7])
				upper = float(parameters[par_keys[par]][8])

				if parameters[par_keys[par]][9] is None and opts['mut_rate'] >= custom.random.random():
					if parameters[par_keys[par]][6] == 'uniform':
						population[par_keys[par], ind] = custom.random.uniform(lower, upper)
					if parameters[par_keys[par]][6] == 'loguniform':
						population[par_keys[par], ind] = numpy.exp(custom.random.uniform(numpy.log(lower), numpy.log(upper)))

				elif parameters[par_keys[par]][9] is not None and float(parameters[par_keys[par]][9]) >= custom.random.random():
					if parameters[par_keys[par]][6] == 'uniform':
						population[par_keys[par], ind] = custom.random.uniform(lower, upper)
					if parameters[par_keys[par]][6] == 'loguniform':
						population[par_keys[par], ind] = numpy.exp(custom.random.uniform(numpy.log(lower), numpy.log(upper)))

	return population

def clean():
	filelist = []
	fileregex = [
		'log*txt',    # log file
		'*.bngl',     # bng2 simulation files
		'*.xml',      # nfsim simulation files. Produced by bng2
		'*.rnf',      # nfsim configuration files
		'*.gdat',     # bng2 simulation outputs (SSA and others)
		'*.cdat',     # bng2 simulation outputs (CVODE)
		'*.species',  # bng2 network generation outputs
		'*.kappa',    # kasim original model and piskas simulations files.
		'model*.sh',  # kasim configuration files
		'model*.txt', # kasim, piskas simulation outputs. Also calculation error outputs
		opts['bin_file'],          # kasim compiled model (no longer done)
		opts['outfile'] + '*.txt', # summary per iteration
	]

	for regex in fileregex:
		filelist.append(glob.glob(regex))
	filelist = [ item for sublist in filelist for item in sublist ]

	for filename in filelist:
		if filename not in [ opts['model'], opts['crit_vals'] ]:
			os.remove(filename)
	#cmd = 'find . -maxdepth 1 -type f ( -name log*txt -o -name *.bngl -o -name *.xml -o -name *.gdat -o -name *.species -o -name *.rnf -o -name model*.sh -o -name model*.out.txt -o -name model*.txt -o -name {:s}*.txt -o -name {:s} ) ! -name {:s} ! -name {:s} -exec rm {{}} +'.format(opts['outfile'], opts['bin_file'], opts['model'], opts['crit_vals'])
	#cmd = re.findall(r'(?:[^\s,"]|"+(?:=|\\.|[^"])*"+)+', cmd)
	#out, err = subprocess.Popen(cmd, shell = False, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()

	return 0

def backup():
	folders = {
		'results' : opts['results'] + '_' + opts['systime'],
		'ranking' : opts['results'] + '_' + opts['systime'] + '/' + opts['ranking'],
		'parsets' : opts['results'] + '_' + opts['systime'] + '/' + opts['parsets'],
		'rawdata' : opts['results'] + '_' + opts['systime'] + '/' + opts['rawdata'],
		'fitness' : opts['results'] + '_' + opts['systime'] + '/' + opts['fitness'],
	}

	# make backup folders
	for folder in folders.values():
		os.mkdir(folder)
	#cmd = 'mkdir -p {:s} {:s} {:s} {:s}'.format(folder_ranking, folder_parsets, folder_rawdata, folder_fitness)
	#cmd = re.findall(r'(?:[^\s,"]|"+(?:=|\\.|[^"])*"+)+', cmd)
	#out, err = subprocess.Popen(cmd, shell = False, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()

	# archive ranking files
	filelist = glob.glob('{:s}*.txt'.format(opts['outfile']))
	for filename in filelist:
		shutil.move(filename, folders['ranking'])
	#cmd = 'find . -maxdepth 1 -type f -name {:s}*.txt -exec mv -t {:s} {{}} +'.format(opts['outfile'], folder_ranking)
	#cmd = re.findall(r'(?:[^\s,"]|"+(?:=|\\.|[^"])*"+)+', cmd)
	#out, err = subprocess.Popen(cmd, shell = False, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()

	# archive simulation outputs
	filelist = glob.glob('model_*.out.txt')
	for filename in filelist:
		shutil.move(filename, folders['rawdata'])
	#cmd = 'find . -maxdepth 1 -type f -name model_*.out.txt -exec mv -t {:s} {{}} +'.format(folder_rawdata)
	#cmd = re.findall(r'(?:[^\s,"]|"+(?:=|\\.|[^"])*"+)+', cmd)
	#out, err = subprocess.Popen(cmd, shell = False, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()

	# archive simulated models
	filelist = glob.glob('model_*.sh')
	for filename in filelist:
		shutil.move(filename, folders['parsets'])
	#cmd = 'find . -maxdepth 1 -type f -name model_*.sh -exec mv -t {:s} {{}} +'.format(folder_parsets)
	#cmd = re.findall(r'(?:[^\s,"]|"+(?:=|\\.|[^"])*"+)+', cmd)
	#out, err = subprocess.Popen(cmd, shell = False, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()

	# archive goodness of fit outputs
	filelist = glob.glob('model_*.txt')
	for filename in filelist:
		shutil.move(filename, folders['fitness'])
	#cmd = 'find . -maxdepth 1 -type f -name model_*.txt -exec mv -t {:s} {{}} +'.format(folder_fitness)
	#cmd = re.findall(r'(?:[^\s,"]|"+(?:=|\\.|[^"])*"+)+', cmd)
	#out, err = subprocess.Popen(cmd, shell = False, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()

	# archive a log file
	log_file = 'log_{:s}.txt'.format(opts['systime'])
	with open(log_file, 'w') as file:
		file.write('# Output of {:s} {:s}\n'.format(opts['python'], subprocess.list2cmdline(sys.argv[0:])))
		#for arg in vars(args): # under test
			#key, value = getattr(args, arg)
			#file.write('{:s}={:s}'.format(key, value))

	shutil.move(log_file, folders['results'])
	shutil.move(opts['model'], folders['results'])
	#cmd = 'cp {:s} {:s} {:s}'.format(opts['model'], log_file, folder_results)
	#cmd = re.findall(r'(?:[^\s,"]|"+(?:=|\\.|[^"])*"+)+', cmd)
	#out, err = subprocess.Popen(cmd, shell = False, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()

	return 0

if __name__ == '__main__':
	sbatch_error = b'sbatch: error: slurm_receive_msg: Socket timed out on send/recv operation\n' \
		b'sbatch: error: Batch job submission failed: Socket timed out on send/recv operation'
	squeue_error = b'squeue: error: slurm_receive_msg: Socket timed out on send/recv operation\n' \
		b'slurm_load_jobs error: Socket timed out on send/recv operation'
	#sbatch_error = b'sbatch: error: Slurm temporarily unable to accept job, sleeping and retrying.'
	#sbatch_error = b'sbatch: error: Batch job submission failed: Resource temporarily unavailable'

	# general options
	args = argsparser()
	opts = ga_opts()
	seed = custom.random.seed(opts['rng_seed'])

	# perform safe checks prior to any calculation
	safe_checks()

	# clean the working directory
	clean()

	# read model configuration
	parameters = configurate()
	# generate first population
	population = populate()
	# generate simulation package
	#population = make_bin() # -mode batch behavior changed from the beta version.

	# main Genetic Algorithm
	for iter in range(opts['num_iter']):
		population = simulate()
		population = evaluate()
		population = ranking()
		population = mutate()

	# move and organize results
	backup()
