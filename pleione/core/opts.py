# -*- coding: utf-8 -*-

'''
Project "Genetic Algorithm for Rule-Based Models", Rodrigo Santibáñez, 2017 @ Dlab, FCV (rsantibanez@dlab.cl)
A Genetic Algorithm inspired by Alberto Martin's Genetic Algorithm, 2016 @ Dlab, FCV (ajmm@dlab.cl)
Citation:
'''

__author__  = 'Rodrigo Santibáñez'
__license__ = 'gpl-3.0'

import numpy, os, time

def ga_opts(args, **kwargs):
	dct = {
		# user defined options
		'seed'      : args.seed,
		'model'     : args.model[0],
		'error'     : args.error,
		'data'      : args.data,
		'slurm'     : args.slurm,
		'others'    : args.sbatch,
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
		'outfile'   : args.output,
		'results'   : args.results,
		'parsets'   : args.parsets,
		'rawdata'   : args.rawdata,
		'fitness'   : args.fitness,
		'ranking'   : args.ranking,
		'python'    : os.path.expanduser(args.python),
		# non-user defined options
		'home'      : os.getcwd(),
		'null'      : '/dev/null',
		'max_error' : numpy.nan,
		'systime'   : str(time.time()).split('.')[0],
		# DEPRECATED just for publication
		'legacy'    : args.legacy,
		'dev'       : args.dev,
		}

	# specific options per simulator
	if kwargs['simulator'] == 'BNG2':
		dct['bng2'] = os.path.expanduser(args.bng2)
	if not kwargs['simulator'] == 'BNG2':
		dct['final'] = args.final[0]
		dct['steps'] = args.steps[0]
	if kwargs['simulator'] == 'KaSim v4':
		dct['kasim'] = os.path.expanduser(args.kasim)
		dct['syntax'] = args.syntax
	if kwargs['simulator'] == 'KaSim v4-beta':
		dct['kasim'] = os.path.expanduser(args.kasim)
		dct['binary'] = args.binary
		dct['bin_file'] = args.model[0].split('.')[0] + '.bin'
	if kwargs['simulator'] == 'NFsim':
		dct['nfsim'] = os.path.expanduser(args.nfsim)
		dct['equil'] = args.equil
	if kwargs['simulator'] == 'PISKaS':
		dct['piskas'] = os.path.expanduser(args.piskas)
		dct['sync'] = args.sync
	dct['simulator'] = kwargs['simulator']

	return dct
