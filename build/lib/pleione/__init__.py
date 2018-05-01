#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''
Project Genetic Algorithm, Rodrigo Santibáñez, 2017-2018 @ Dlab, FCV, Chile (rsantibanez@dlab.cl)
A Python implementation of Alberto Martin's Genetic Algorithm, 2016 @ Dlab, FCV (proteinomano@gmail.cl)
To be used with NFsim. Please refer to other subprojects for other stochastic simulators support.
Citation:
'''

__author__  = 'Rodrigo Santibáñez <rsantibanez@dlab.cl>, Computational Biology Lab, Fundacion Ciencia y Vida, Chile'
__license__ = 'gpl-3.0'
__version__ = 'final-nfsim-v1.12.1'

import argparse, glob, multiprocessing, os, re, shlex, subprocess, time
import pandas, numpy

def parallelize(cmd):
        print(cmd)
        proc = subprocess.Popen(shlex.split(cmd), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        out, err = proc.communicate()
        proc.wait()

        return 0

