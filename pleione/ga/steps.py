# -*- coding: utf-8 -*-

'''
Project "Genetic Algorithm for Rule-Based Models", Rodrigo Santibáñez, 2017 @ Dlab, FCV (rsantibanez@dlab.cl)
A Genetic Algorithm inspired by Alberto Martin's Genetic Algorithm, 2016 @ Dlab, FCV (ajmm@dlab.cl)
Citation:
'''

__author__  = 'Rodrigo Santibáñez'
__license__ = 'gpl-3.0'

import glob, multiprocessing, numpy, os, pandas, re, subprocess, sys, time
from pleione.core.utils import parallelize, write_models
from pleione.core.random import random

def simulate(opts):
	# slice dictionary
	parameters = opts['parameters']
	population = opts['population']

	job_desc = {
		'nodes'     : 1,
		'ntasks'    : 1,
		'ncpus'     : 1,
		'null'      : opts['null'],
		'partition' : opts['slurm'],
		'others'    : opts['others'],
		'job_name'  : 'child_{:s}'.format(opts['systime']),
		'stdout'    : 'stdout_{:s}.txt'.format(opts['systime']),
		'stderr'    : 'stderr_{:s}.txt'.format(opts['systime']),
		}

	# write models
	write_models(opts)

	# submit simulations to the queue
	squeue = []
	model_string = '{:s}.{:0' + str(len(str(opts['num_sims']))) + 'd}.out.txt'
	for ind in range(opts['pop_size']):
		for sim in range(opts['num_sims']):
			model = population['model', ind]
			output = model_string.format(model, sim)

			if not os.path.exists(output):
				job_desc['exec_kasim'] = '{:s} -i {:s}.kappa -l {:s} -p {:s} -o {:s} -syntax {:s} --no-log'.format( \
					opts['kasim'], model, opts['final'], opts['steps'], output, opts['syntax'])

				# use SLURM Workload Manager
				if opts['slurm'] is not None:
					cmd = os.path.expanduser('sbatch --no-requeue -p {partition} -N {nodes} -c {ncpus} -n {ntasks} -o {null} -e {null} -J {job_name} {others} \
						--wrap ""{exec_kasim}""'.format(**job_desc))
					cmd = re.findall(r'(?:[^\s,"]|"+(?:=|\\.|[^"])*"+)+', cmd)
					out, err = subprocess.Popen(cmd, shell = False, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
					while err == sbatch_error:
						out, err = subprocess.Popen(cmd, shell = False, stdout = subprocess.PIPE, stderr = subprocess.PIPE).communicate()
					squeue.append(out.decode('utf-8')[20:-1])

				# to use with multiprocessing.Pool
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

def evaluate(opts):
	# slice dictionary
	population = opts['population']

	job_desc = {
		'nodes'     : 1,
		'ntasks'    : 1,
		'ncpus'     : 1,
		'null'      : opts['null'],
		'partition' : opts['slurm'],
		'job_name'  : 'child_{:s}'.format(opts['systime']),
		'stdout'    : 'stdout_{:s}.txt'.format(opts['systime']),
		'stderr'    : 'stderr_{:s}.txt'.format(opts['systime']),
		'doerror'   : '{:s} -m pleione.kasim-doerror '.format(opts['python']),
		}

	if set(opts['error']).issuperset(set(['MWUT'])) or opts['dev']:
		job_desc['doerror'] = job_desc['doerror'] + '--crit {:s} '.format(opts['crit_vals'])
		job_desc['deverror'] = job_desc['deverror'] + '--crit {:s} '.format(opts['crit_vals'])

	# submit error calculations to the queue
	squeue = []

	for ind in range(opts['pop_size']):
		model = population['model', ind]

		# subscript cannot expand asterisk wildcard (for security reasons)
		# also, Popen interprets it as an literal asterisk, making necesary the interpretation with glob
		# this makes a problem when pleione is a subroutine or the main routine.
		# Attempt to reconcile both situations:
		data = ''
		for value in opts['data']:
			if '*' in value:
				data += ' '.join(glob.glob(value)) + ' '
			else:
				data += value + ' '

		error = ' '.join(opts['error'])
		sims = ' '.join(glob.glob('{:s}.*.out.txt'.format(model)))
		output = '{:s}.txt'.format(model)

		job_desc['calc'] = job_desc['doerror'] + '--data {:s} --sims {:s} --file {:s} --error {:s}'.format(data, sims, output, error)
		if opts['dev']:
			job_desc['calc'] = job_desc['doerror'] + '--data {:s} --sims {:s} --file {:s} --do_all True'.format(data, sims, output)

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

def ranking(opts):
	# slice dictionary
	parameters = opts['parameters']
	population = opts['population']

	for ind in range(opts['pop_size']):
		with open('{:s}.txt'.format(population['model', ind]), 'r') as file:
			tmp = pandas.read_csv(file, delimiter = '\t', header = None)
			data = tmp.set_index(0, drop = False).rename_axis(None, axis = 0).drop(0, axis = 1).rename(columns = {1: 'value'})

		if opts['dev']:
			fitfunc = list(data.index)
		else:
			fitfunc = opts['error']

		# and store the error in the population dictionary
		for name in range(len(fitfunc)):
			population[fitfunc[name], ind] = data.loc[fitfunc[name], 'value']

	# now that everything is stored in the population dict, we proceed to rank models by the selected error function(s)
	jobs = {}
	rank = {}
	fitfunc = opts['error']

	for name in range(len(fitfunc)):
		for ind in range(opts['pop_size']):
			jobs[population['model', ind]] = population[fitfunc[name], ind]
		rank[fitfunc[name]] = sorted(jobs, key = jobs.get, reverse = False)

	for name in range(len(fitfunc)):
		for ind in range(opts['pop_size']):
			jobs[population['model', ind]] += { key : value for value, key in enumerate(rank[fitfunc[name]]) }[population['model', ind]]

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

	if opts['dev']:
		fitfunc = sorted(list(data.index))

	par_string = '{:.' + opts['par_fmt'] + '}\t'
	iter_string = '{:s}_{:0' + str(len(str(opts['num_iter']))) + 'd}.txt'
	with open(iter_string.format(opts['outfile'], opts['iter']), 'w') as file:
		file.write('# Output of {:s} {:s}\n'.format(opts['python'], subprocess.list2cmdline(sys.argv[0:])))
		file.write('Elapsed time: {:.0f} seconds\n'.format(time.time() - float(opts['systime'])))
		file.write('iteration: {:03d}\t'.format(opts['iter']))
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
				if fitfunc[i] not in ['DUT', 'MWUT', 'WMWET', 'TOST']:
					file.write(par_string.format(float(population[fitfunc[i], ind])))
				else:
					file.write('{:.0f}\t'.format(float(population[fitfunc[i], ind])))
			for key in range(len(par_keys)):
				file.write(par_string.format(float(population[par_keys[key], ind])))
			file.write('\n')

	return population

def mutate(opts):
	# slice dictionary
	parameters = opts['parameters']
	population = opts['population']

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
			if not opts['legacy']:
				# choose two random individuals from the best population (reindexed from 0 to 'pop_best' size)
				n1 = numpy.random.choice(range(top), p = prob[0:top])
				n2 = numpy.random.choice(range(top), p = prob[0:top])
			else:
				n1 = random.choice(range(top))
				n2 = random.choice(range(top))

			if opts['self_rec'] == False and not opts['pop_size'] == 1:
				while n2 == n1:
					if opts['legacy']:
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
				if opts['mut_swap'] >= random.random(opts['legacy']):
					population[par_keys[par], ind] = best_population[par_keys[par], n2]
					population[par_keys[par], ind + 1] = best_population[par_keys[par], n1]

		elif opts['xpoints'] == 'single':
			point = random.uniform(0, len(par_keys))
			for par in range(len(par_keys)):
				# create children and do not swap parameters!
				if par <= point:
					population[par_keys[par], ind] = best_population[par_keys[par], n1]
					population[par_keys[par], ind + 1] = best_population[par_keys[par], n2]
				else:
					population[par_keys[par], ind] = best_population[par_keys[par], n2]
					population[par_keys[par], ind + 1] = best_population[par_keys[par], n1]

		# include the model id
		model_string = 'model_{:0' + str(len(str(opts['num_iter']))) + 'd}' + '_{:0' + str(len(str(opts['pop_size']))) + 'd}'
		population['model', ind] = model_string.format(opts['iter'] + 1, ind)
		population['model', ind + 1] = model_string.format(opts['iter'] + 1, ind + 1)

		# include the error in the population dictionary
		population['error', ind] = opts['max_error']
		population['error', ind + 1] = opts['max_error']

	# mutate parameter values
	for ind in range(opts['pop_best'], opts['pop_size']):
		for par in range(len(par_keys)):
			if parameters[par_keys[par]][6] == 'factor':
				if float(parameters[par_keys[par]][7]) >= random.random(opts['legacy']):
					population[par_keys[par], ind] = population[par_keys[par], ind] * \
						random.uniform(opts['legacy'], 1.0 - float(parameters[par_keys[par]][8]), 1.0 + float(parameters[par_keys[par]][8]))

			elif parameters[par_keys[par]][6] == 'uniform' or parameters[par_keys[par]][6] == 'loguniform':
				lower = float(parameters[par_keys[par]][7])
				upper = float(parameters[par_keys[par]][8])

				if parameters[par_keys[par]][9] is None and opts['mut_rate'] >= random.random(opts['legacy']):
					if parameters[par_keys[par]][6] == 'uniform':
						population[par_keys[par], ind] = random.uniform(opts['legacy'], lower, upper)
					if parameters[par_keys[par]][6] == 'loguniform':
						population[par_keys[par], ind] = numpy.exp(random.uniform(opts['legacy'], numpy.log(lower), numpy.log(upper)))

				elif parameters[par_keys[par]][9] is not None and float(parameters[par_keys[par]][9]) >= random.random():
					if parameters[par_keys[par]][6] == 'uniform':
						population[par_keys[par], ind] = random.uniform(opts['legacy'], lower, upper)
					if parameters[par_keys[par]][6] == 'loguniform':
						population[par_keys[par], ind] = numpy.exp(random.uniform(opts['legacy'], numpy.log(lower), numpy.log(upper)))

	return population
