Parameterization with NFsim
===========================

1. **Prepare the model**

   Pleione finds which variables will be calibrated using
   the symbol ``#`` (number sign, hash or pound sign) followed by:

	* An initial distribution type: ``uniform``, ``loguniform``, ``lognormal``
	* An initial search space: ``[min max]`` or ``[mean standard_deviation]``
	  in the case if ``lognormal`` was selected.
	* A type of mutation: ``uniform`` or ``loguniform`` to use a new search
	  space; or ``factor`` to perform a local mutation search
	* A search space for mutated parameters: ``[min max]`` or
	  ``[probability fold_change]``
	* An optional mutation rate per parameter. Without it, a global mutation
	  rate is used.

   For instace:

.. code-block:: bash

	  KD1__FREE__        1.000000e+00 # loguniform[0.01 100] factor[0.2 0.1]
	  km1__FREE__        1.000000e+00 # loguniform[0.01 100] factor[0.2 0.1]
	  K2RT__FREE__       1.000000e+00 # loguniform[0.01 100] factor[0.2 0.1]
	  km2__FREE__        1.000000e+00 # loguniform[0.01 100] factor[0.2 0.1]
	  kphos__FREE__      1.000000e+00 # loguniform[0.01 100] factor[0.2 0.1]
	  kdephos__FREE__    1.000000e+00 # loguniform[0.01 100] factor[0.2 0.1]

.. note::
	**Factor mutation:** This type of mutation strategy comes from BioNetFit and
	selects a random value from the range ``0.9 * old_value, 1.1 * old_value``
	if the declared value is ``0.1`` with probability ``0.2``.

2. **Prepare the data files**

   NFsim produce simulations files with the following format. Please prepare
   data files with the same format to avoid incompatibilities.

.. code-block:: bash

	time, RLbonds, pR
	0.00000000E+00, 0.00000000E+00, 3.55300000E+02
	1.00000000E+01, 1.14072000E+02, 3.56440000E+02
	2.00000000E+01, 1.39183800E+02, 3.49960000E+02
	3.00000000E+01, 1.49153400E+02, 3.43980000E+02
	4.00000000E+01, 1.56868400E+02, 3.42600000E+02
	5.00000000E+01, 1.56788000E+02, 3.35620000E+02
	6.00000000E+01, 1.63666800E+02, 3.37480000E+02

2. **Prepare a sbatch configuration file**

   Use the following code as template to make a shell script and queue it with
   sbatch. Note that the ``export`` statement is inside the code to tell SLURM
   to add the path and ensure proper execution when pleione was cloned with
   git. Also, ``python3`` redirects to either the system installed executable
   (with pandas installed either as admin or user) or redirects to the user
   compiled executable if an alias exists for it.

.. code-block:: bash

	#!/bin/sh

	#SBATCH --no-requeue
	#SBATCH --partition=cpu

	#SBATCH --nodes=1
	#SBATCH --ntasks=1
	#SBATCH --cpus-per-task=1

	#SBATCH --job-name=pleione-nfsim
	#SBATCH --output=stdout.txt
	#SBATCH --error=stderr.txt

	export PYTHONPATH="$PYTHONPATH:$HOME/opt/git-glucksfall-pleione-master"

	MODEL=pysbmodel-example6-nfsim.bngl # the model should have the .bngl extension
	FINAL=60
	STEPS=6

	PARTITION=$SLURM_JOB_PARTITION
	DATA=../exp-data/nfsim/data-*.txt

	NUM_ITER=100
	NUM_SIMS=10
	POP_SIZE=100
	POP_BEST=0

	SWAP=0.5
	RATE=0.5
	ERROR="MWUT"
	UTABLE=./ucrit.txt

	python3 -m pleione.nfsim --model=$MODEL --final=$FINAL --steps=$STEPS \
	--iter=$NUM_ITER --pops=$POP_SIZE --sims=$NUM_SIMS --best=$POP_BEST \
	--data=$DATA --rate=$RATE --swap=$SWAP --error=$ERROR --crit=$UTABLE \
	--slurm=$PARTITION

.. note::
	**sbatch or python multiproccesing?** To execute Pleione outside a SLURM
	queue, simple execute the shell script with ``sh``, ``bash`` or any shell
	interpreter without the ``slurm`` option. Be aware that, if SLURM is
	running in the same machine, Pleione subprocess would impact negatively in
	other user's threads, and viceversa, since a cpu core could execute
	concurrently two threads.

.. note::
	**Need help?** type ``python3 -m pleione.nfsim --help`` to find out the
	available command options.

.. refs
.. _KaSim: https://github.com/Kappa-Dev/KaSim
.. _NFsim: https://github.com/RuleWorld/nfsim
.. _BioNetGen2: https://github.com/RuleWorld/bionetgen
.. _PISKaS: https://github.com/DLab/PISKaS
.. _BioNetFit: https://github.com/RuleWorld/BioNetFit
.. _SLURM: https://slurm.schedmd.com/

.. _Kappa: https://www.kappalanguage.org/
.. _BioNetGen: http://www.csb.pitt.edu/Faculty/Faeder/?page_id=409
.. _pandas: https://pandas.pydata.org/
