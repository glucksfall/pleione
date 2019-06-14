.. _Fitneess_Functions:

Objective Functions
===================

Common to all parameterization methods, there are 9 algebraic objective
functions and one statistical function already implemented in the code.
Moreover, the code sort the models by their rank and therefore, ranks can be
added and sorted again, making the possibility to use a Multiple Objective
Genetic Algorithm.

.. toctree::
	:maxdepth: 3

	functions/algebraics
	functions/statistical
	functions/multiple
	functions/addFitness

.. note::
	**Need a different Objective Function?**
	The code that calculates the error is separated from the main Genetic
	Algorithm. This make useful to encode other Objective Functions if the
	already implemented does not apply to your necessities. You could contact us
	to add your function to the pleione package.
