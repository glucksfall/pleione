Statistical Objective Functions
===============================

We implemented the Mann-Whitney U-test (MWUT) to calculate the error between
experimental data and simulations. The U-test is a non-parametric statistical
test that, within a confidence level, determine if a stochastic repeated
measurements is identical or not to another repeated measurements.

1. We count how many times experimental data (:math:`exp_i`) are larger than
simulated values (:math:`sim_j`):

	| **for** :math:`i \mathrm{\ in\ } \mathrm{range} ( \mathrm{len}(exp) )`:
	|   **for** :math:`j \mathrm{\ in\ } \mathrm{range} ( \mathrm{len}(sim) )`:
	|      **if** :math:`exp_{i} > sim_{j}`:
	|         :math:`U_{exp} \gets U_{exp} + 1.0`
	|      **else if** :math:`exp_{i} < sim_{j}`:
	|         :math:`U_{sim} \gets U_{sim} + 1.0`
	|      **else**:
	|         :math:`U_{exp} \gets U_{exp} + 0.5`
	|         :math:`U_{sim} \gets U_{sim} + 0.5`

2. We determine if :math:`U_{exp}` is statistically significant:

	| **for** :math:`i \mathrm{\ in\ } \mathrm{range} ( \mathrm{len}(exp) )`:
	|   **if** :math:`\mathrm{len}(exp) \times \mathrm{len}(sim) - \mathrm{min}(U_{exp}, U_{sim}) \leq U_{\mathrm{critic}}`:
	|      :math:`\mathrm{\textit{null}\ hypothesis,\ }H_{0}\mathrm{,\ is\ rejected}`
	|      :math:`U_{\mathrm{model}} \gets U_{\mathrm{model}} + 1.0`

.. note::
	The U-test is the only fitness function that has known limits: For a
	*perfect* model, the U-test is zero. A complete wrong model will have a
	:math:`U_{model}` equal to the number of Observables times the number of
	experimental time points. For instance, the example model we use to compare
	with BioNetFit has 2 Observables and 7 experimental time points, then a max
	:math:`U_{model}` equal to 14.
