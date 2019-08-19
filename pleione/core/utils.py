# -*- coding: utf-8 -*-

'''
Project "Genetic Algorithm for Rule-Based Models", Rodrigo Santibáñez, 2017 @ Dlab, FCV (rsantibanez@dlab.cl)
A Genetic Algorithm inspired by Alberto Martin's Genetic Algorithm, 2016 @ Dlab, FCV (ajmm@dlab.cl)
Citation:
'''

__author__  = 'Rodrigo Santibáñez'
__license__ = 'gpl-3.0'

import glob, os, shutil, subprocess, sys

# Catch errors from SLURM
"""
sbatch_error = b'sbatch: error: slurm_receive_msg: Socket timed out on send/recv operation\n' \
	b'sbatch: error: Batch job submission failed: Socket timed out on send/recv operation'
squeue_error = b'squeue: error: slurm_receive_msg: Socket timed out on send/recv operation\n' \
	b'slurm_load_jobs error: Socket timed out on send/recv operation'
sbatch_error = b'sbatch: error: Slurm temporarily unable to accept job, sleeping and retrying.'
sbatch_error = b'sbatch: error: Batch job submission failed: Resource temporarily unavailable'
"""

def parallelize(cmd):
	proc = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
	out, err = proc.communicate()
	proc.wait()
	return 0

def write_models(opts):
	# slice dictionary
	parameters = opts['parameters']
	population = opts['population']

	if opts['simulator'] == 'KaSim v4':
		# generate a new kappa file per model
		par_keys = list(parameters.keys())
		par_string = '%var: \'{:s}\' {:.' + opts['par_fmt'] + '}\n'
		for ind in range(opts['pop_size']):
			model = population['model', ind]

			if not os.path.exists(model + '.kappa'):
				with open(model + '.kappa', 'w') as file:
					for line in range(len(par_keys)):
						if parameters[line][0] == 'par':
							file.write(par_string.format(parameters[line][1], population[line, ind]))
						else:
							file.write(parameters[line])

	return 0

def clean(opts):
	filelist = []
	fileregex = [
		'log*txt',     # log file
		'*bngl',       # bng2 simulation files
		'*xml',        # nfsim simulation files. Produced by bng2
		'*rnf',        # nfsim configuration files
		'*gdat',       # bng2 simulation outputs (SSA and others)
		'*cdat',       # bng2 simulation outputs (CVODE)
		'*species',    # bng2 network generation outputs
		'*kappa',      # kasim original model and piskas simulations files.
		'model_*.txt', # kasim, piskas simulation outputs.
		               # Also calculation error outputs
		opts['outfile'] + '*.txt', # summaries per iteration

		#DEPRECATED: KaSim v4-beta
		#'model_*.sh',     # kasim configuration files.
		#opts['bin_file'], # kasim compiled model.
	]

	for regex in fileregex:
		filelist.append(glob.glob(regex))
	filelist = [ item for sublist in filelist for item in sublist ]

	for filename in filelist:
		if filename not in [ opts['model'], opts['crit_vals'] ]:
			os.remove(filename)

	return 0

def backup(opts):
	results = opts['results'] + '_' + opts['systime']
	folders = {
		'ranking' : results + '/' + opts['ranking'],
		'parsets' : results + '/' + opts['parsets'],
		'rawdata' : results + '/' + opts['rawdata'],
		'fitness' : results + '/' + opts['fitness'],
	}

	# make backup folders
	os.mkdir(results)
	for folder in folders.values():
		os.mkdir(folder)

	# archive ranking files
	filelist = glob.glob('{:s}*.txt'.format(opts['outfile']))
	for filename in filelist:
		shutil.move(filename, folders['ranking'])

	# archive simulation outputs
	if opts['simulator'] == 'KaSim v4':
		filelist = glob.glob('model_*.out.txt')

	for filename in filelist:
		shutil.move(filename, folders['rawdata'])

	# archive models
	if opts['simulator'] == 'KaSim v4':
		filelist = glob.glob('model_*.kappa')

	for filename in filelist:
		shutil.move(filename, folders['parsets'])

	# archive goodness of fit outputs
	filelist = glob.glob('model_*.txt')
	for filename in filelist:
		shutil.move(filename, folders['fitness'])

	# archive a log file
	log_file = 'log_{:s}.txt'.format(opts['systime'])
	with open(log_file, 'w') as file:
		file.write('# Output of {:s} {:s}\n'.format(opts['python'], subprocess.list2cmdline(sys.argv[0:])))
	shutil.move(log_file, results)
	shutil.copy2(opts['model'], results)

	clean(opts)

	return 0
