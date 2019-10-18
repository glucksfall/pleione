# -*- coding: utf-8 -*-

'''
Project "Genetic Algorithm for Systems Biology", Rodrigo Santibáñez, 2019 @ NBL, Universidad Mayor
An extension of Pleione v1.5: Pleione: A tool for statistical and multi-objective calibration of Rule-based models. Scientific Reports (2019)
To be used with Tellurium v2.1.5. Please refer to other subprojects for other stochastic simulators support
'''

__author__  = 'Rodrigo Santibáñez'
__license__ = 'gpl-3.0'
__software__ = 'tellurium-v2.1.5'

import os, sys
import pandas, tellurium

rr = tellurium.loada(sys.argv[1])
rr.integrator = 'rk4'
results = rr.simulate(0, float(sys.argv[2]), int(sys.argv[3]) + 1)

with open(sys.argv[4], 'w') as outfile:
	pandas.DataFrame(data = results, columns = results.colnames).fillna(0).to_csv(outfile, sep = '\t', index = False)
