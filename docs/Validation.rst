Model Validation
================

Pleione's parameter calibration scripts call an external script to calculate
fitness to experimental data. You could use one of the following script to
calculate the fitness of your parameterized model against an independent
experimental data set:

.. code-block:: bash

	python3 -m pleione.bng2-doerror --data foo --sims bar \
	--file output.txt --error SDA

*OR*

.. code-block:: bash

	python3 -m pleione.kasim-doerror --data foo --sims bar \
	--file output.txt --error SDA

*OR*

.. code-block:: bash

	python3 -m pleione.nfsim-doerror --data foo --sims bar \
	--file output.txt --error SDA

*OR*

.. code-block:: bash

	python3 -m pleione.piskas-doerror --data foo --sims bar \
	--file output.txt --error SDA

.. note::
	**Fitness Function**

	Pleione currently support 9 algebraics and 3 statistical tests as fit functions.
	To calculate more than one function, include a list such as ``SDA SSQ CHISQ``.
	In doing so, this will calculate the contribution of both o more metrics to
	the overall error and aid to validate of discard a model calibration. More information
	in :ref:`Fitneess_Functions`

.. note::
	**All fitness functions**

	Use the ``--do_all True`` argument to calculate all fitness functions included in Pleione.
	You should provide a table with critical values (see example folder for the one-tail table)

.. note::
	**(non-)Rejection matrices**

	The ``--report True`` argument will print to the console relevant calculations for the
	statistical tests.

.. note::
	**Need Help?**
	Type ``python3 -m pleione.$STOCH_ENGINE-doerror --help`` where
	``$STOCH_ENGINE`` can be the currently supported stochastic engines: BNG2,
	NFsim, KaSim and PISKaS (all in lower cases, for instance ``nfsim-doerror``)
