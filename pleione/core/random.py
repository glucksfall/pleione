# -*- coding: utf-8 -*-

'''
Project "Genetic Algorithm for Rule-Based Models", Rodrigo Santibáñez, 2017 @ Dlab, FCV (rsantibanez@dlab.cl)
A Genetic Algorithm inspired by Alberto Martin's Genetic Algorithm, 2016 @ Dlab, FCV (ajmm@dlab.cl)
Citation:
'''

__author__  = 'Rodrigo Santibáñez'
__license__ = 'gpl-3.0'

import random, numpy

class random:
	def seed(n, legacy):
		if legacy:
			random.seed(n)
		else:
			numpy.random.seed(n)

	def random(legacy):
		if legacy:
			return random.random()
		else:
			return numpy.random.random()

	def uniform(legacy, lower, upper):
		if legacy:
			return random.uniform(lower, upper)
		else:
			return numpy.random.uniform(lower, upper, None)

	def lognormal(legacy, lower, upper):
		if legacy:
			return random.lognormvariate(lower, upper)
		else:
			return numpy.random.lognormal(lower, upper, None)
