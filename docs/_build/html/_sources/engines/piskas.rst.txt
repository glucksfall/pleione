Parameterization with PISKaS
============================

1. **Prepare the model**

   Pleione parmeterization methods find which variables will be calibrated using
   the symbol ``#`` (number sign, hash or pound sign) followed by:

	* An initial distribution type: ``uniform``, ``loguniform``, ``lognormal``
	* An initial search space: ``[min max]`` or ``[mean standard_deviation]``
	  in the case ``lognormal`` was selected.
	* A type of mutation: ``uniform`` or ``loguniform`` to use a new search
	  space; or ``factor`` to perform a local mutation search
	* A search space for mutated parameters: ``[min max]`` or
	  ``[probability fold_change]``
	* An optional mutation rate per parameter. Without it, a global mutation
	  rate is used.

   For instace:

.. code-block:: bash

	%var: 'KD1__FREE__' 1.000000e+00 # loguniform[0.01 100] factor[0.2 0.1]
	%var: 'km1__FREE__' 1.000000e+00 # loguniform[0.01 100] factor[0.2 0.1]
	%var: 'K2RT__FREE__' 1.000000e+00 # loguniform[0.01 100] factor[0.2 0.1]
	%var: 'km2__FREE__' 1.000000e+00 # loguniform[0.01 100] factor[0.2 0.1]
	%var: 'kphos__FREE__' 1.000000e+00 # loguniform[0.01 100] factor[0.2 0.1]
	%var: 'kdephos__FREE__' 1.000000e+00 # loguniform[0.01 100] factor[0.2 0.1]

.. note::
	**Factor mutation:** This type of mutation strategy comes from BioNetFit and
	selects a random value from the range ``0.9 * old_value, 1.1 * old_value``
	if the declared value is ``0.1`` with probability ``0.2``.

2. **Prepare the data files**

   PISKaS produce simulations files with the following format. Please prepare
   data files with the same format, including the initial space in each row.

.. code-block:: bash

	example time 'RLbonds' 'pR'
	 6.000000E+02 0 355.3
	 6.100000E+02 114.072 356.44
	 6.200000E+02 139.1838 349.96
	 6.300000E+02 149.1534 343.98
	 6.400000E+02 156.8684 342.6
	 6.500000E+02 156.788 335.62
	 6.600000E+02 163.6668 337.48

.. note::
	**An extra column name?** PISKaS produces one output for each compartment
	declared in the model. Therefore, adding the name of the compartment as the
	first column allows the code to identify the corresponding experimental data
	with the simulated compartment and apply correctly the fitness function.
	Finally, the model error is the sum of each fitness per compartment.

.. note::
	**About the example model:** The model has three parts: An equilibration of
	600 seconds, then the model is modified to add a quantity of ``L(r)`` agents,
	and then perform the actual simulation for 60 seconds. Despite BNG2 and NFsim,
	PISKaS reports the whole simulation, so to compare effectively, we must offset
	the time of the experimental data by 600.

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

	#SBATCH --job-name=pleione-piskas
	#SBATCH --output=stdout.txt
	#SBATCH --error=stderr.txt

	export PYTHONPATH="$PYTHONPATH:$HOME/opt/git-glucksfall-pleione-master"

	MODEL=pysbmodel-example6-piskas.kappa
	FINAL=660
	STEPS=66

	PARTITION=$SLURM_JOB_PARTITION
	DATA=../exp-data/piskas/data-*.txt

	NUM_ITER=100
	NUM_SIMS=10
	POP_SIZE=100
	POP_BEST=0

	SWAP=0.5
	RATE=0.5
	ERROR="MWUT"
	UTABLE=./ucrit.txt

	python3 -m pleione.piskas --model=$MODEL --final=$FINAL --steps=$STEPS \
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
	**Need help?** type ``python3 -m pleione.piskas --help`` to find out the
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
