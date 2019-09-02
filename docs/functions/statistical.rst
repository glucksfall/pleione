Statistical Objective Functions
===============================

We implemented three equivalence tests to determine if two distribution are
similar in a interval. The interval is defined by default as one standard
deviation of experimental data for the two one-sided t-tests (TOST) and the
Double Mann-Whitney U-test (DUT). In the case of the Wellek's test (WMWET), the
equivalence interval is :math:`\epsilon_1 = 0.3129` and :math:`\epsilon_2 =
0.2661`.

The user can set the ``--factor`` argument to divide the standard deviation by
it, or can set the ``--stdv sims`` argument to use rather the standard
deviation of simulations, or provide custom limits with ``--lower`` and, or
``--upper`` arguments, which point to one file with the same structure as the
experimental data. In the case the user omits ``--lower`` or ``--upper``, the
equivalence interval will be symmetrical.

To calculate TOST, we use the ``ttost_ind`` function from the python
*statsmodels* package. In the case of the Wellek's test, we implemented in
python the ``mawi.R`` script from the EQUIVNONINF package
(https://rdrr.io/cran/EQUIVNONINF/man/mawi.html). And for the Double
Mann-Whitney U-test, we implemented it as two Mann-Whitney U-test as follows:

The U-test is a non-parametric statistical test that, within a confidence level,
determine if a random distribution is different (two-tails) or greater
(one-tail) compared to a second distribution. The Algorithm is valid to compare
distribution of 3 to 20 measurements.

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

	| :math:`U_{\mathrm{model}} = U_{max} = \mathrm{len}(exp) \times \mathrm{len}(sim)`
	| **for** :math:`i \mathrm{\ in\ } \mathrm{range} ( \mathrm{len}(exp) )`:
	|   **for** :math:`j \mathrm{\ in\ } \mathrm{range} ( \mathrm{len}(sim))`:
	|      test :math:`H_0: exp > sim − lower`
	|      **if** :math:`U_{max} - U_{exp} <= U_{critic}` then null hypothesis, :math:`H_0`, is rejected
	|          :math:`U_{lower} = 1.0`
	|      **else**
	|          :math:`U_{lower} = 0.0`
	|      test :math:`H_0: exp < sim + upper`
	|      **if** :math:`U_{max} - U_{sim} <= U_{critic}` then null hypothesis, :math:`H_0`, is rejected
	|          :math:`U_{upper} = 1.0`
	|      **else**
	|          :math:`U_{upper} = 0.0`
	|      :math:`U_{model} = U_{model} - U_{lower} \times U_{upper}`

.. note::
	The iterative statistical tests are fitness functions having known limits: For a
	*perfect* model, the U-test is zero. A complete wrong model will have a
	:math:`U_{model}` equal to the number of Observables times the number of
	experimental time points. For instance, the example model we use to compare
	with BioNetFit has 2 Observables and 7 experimental time points, then a max
	:math:`U_{model}` equal to 14.
