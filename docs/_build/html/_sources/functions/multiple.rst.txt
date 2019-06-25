Multiple Objective Functions
============================

A Multiple Objective Function is build from two or more fitness functions.
Firstly, a fitness is calculated and all models ranked. Then, the next fitness.
Finally, the sum of ranks is use to rank agains the models.

Algoritmically:

.. math::
	rank_1 &= \mathrm{sort\ models\ following\ function\ 1} \\
	&\mathrel{\vdots} \\
	rank_n &= \mathrm{sort\ models\ following\ function\ n} \\
	\\
	rank_{MO} &= \mathrm{sort\ models\ following\ } (rank_1 + \ldots + rank_n)

.. note::
	We currently don't provide weights to rank the models. Be aware that, if you
	use multiple algebraic functions and the statistical fitness function, the
	importance of the statistical function is diluited.
