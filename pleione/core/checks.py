# -*- coding: utf-8 -*-

'''
Project "Genetic Algorithm for Rule-Based Models", Rodrigo Santibáñez, 2017 @ Dlab, FCV (rsantibanez@dlab.cl)
A Genetic Algorithm inspired by Alberto Martin's Genetic Algorithm, 2016 @ Dlab, FCV (ajmm@dlab.cl)
Citation:
'''

__author__  = 'Rodrigo Santibáñez'
__license__ = 'gpl-3.0'

import glob, os, shutil

def safe_checks(opts):
	error_msg = ''
	if shutil.which(opts['python']) is None:
		error_msg += 'python3 (set at {:s}) can\'t be called to perform error calculation.\n' \
			'You could use --python {:s}\n'.format(opts['python'], shutil.which('python3'))

	# check for simulators
	if opts['simulator'] == 'BNG2' and shutil.which(opts['bng2']) is None:
		error_msg += 'BNG2 (set at {:s}) can\'t be called to perform simulations.\n' \
			'Check the path to BNG2.'.format(opts['bng2'])
	if opts['simulator'] == 'KaSim v4' and shutil.which(opts['kasim']) is None:
		error_msg += 'KaSim (set at {:s}) can\'t be called to perform simulations.\n' \
			'Check the path to KaSim.'.format(opts['kasim'])
	if opts['simulator'] == 'NFsim' and shutil.which(opts['nfsim']) is None:
		error_msg += 'NFsim (set at {:s}) can\'t be called to perform simulations.\n' \
			'Check the path to NFsim.'.format(opts['nfsim'])
	if opts['simulator'] == 'PISKaS' and shutil.which(opts['piskas']) is None:
		error_msg += 'PISKaS (set at {:s}) can\'t be called to perform simulations.\n' \
			'Check the path to PISKaS.'.format(opts['piskas'])

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
		for data in opts['data']: # the shell expanded the *
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
