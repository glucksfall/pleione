Parameters estimation
=====================

Pleione's parameterization methods rely on Computational Load Distribution. The
naïve approach is to use the python's *multiprocessing* API and each simulation
distributed within the Pool of available (minus one) cores. This approach would
make pleione's methods compatible with Microsoft Windows and Apple OS X.
However, to take fully advantage of High-Performance Computing architectures,
pleione's methods rely on SLURM --*Simple Linux Utility for Resource
Management*-- (`SLURM`_) to distribute simulations through your infrastructure,
remote infrastructures, and cloud services like Google Compute Engine, Microsoft
Azure, and Amazon Elastic Compute Cloud.

Up to date, pleione's parameterization methods rely on 4 simulations engines:
KaSim and PISKaS simulate *kappa* language models. Unlike KaSim, PISKaS is able
to simulate multiple compartment models distributing the calculation of each
compartment through multiple cores. In the other hand, BioNetGen2 and NFsim
simulate *BioNetGen* language models. Despite KaSim and PISKaS, BioNetGen2 does
not provide a Command-Line Interface to especify simulation parameters and
rather, the simulation parameters (e.g. time to simulation, number of points to
report, ...) must be given inside the model specification. Moreover, you need to
especify the simulation engine to use: Deterministic simulation through *CVODE*,
the Stochastic Simulation Algorithm *SSA*, Exact Hybrid Particle/Population
Algorithm *HPP*, and the Partition-Leap Algorithm *PLA*. Moreover, NFsim could
be used by BioNetGen2 to simulate models or called externally after creating the
model xml especification with BioNetGen2 --xml option.

Because the software requirements and differences, we provide specific
documentation to all of them rather than provide common guidelines and then
stating the differences.

Parameterization of kappa-language Rule-Based Models

.. toctree::
	:maxdepth: 3

	engines/kasim
	engines/piskas

Parameterization of BioNetGen language Rule-Based Models

.. toctree::
	:maxdepth: 3

	engines/bng2
	engines/nfsim

Common to all parameterization methods, there are 9 algebraic objective
functions and one statistical function already implemented in the code.
Moreover, the code sort the models by their rank and therefore, ranks can be
added and sorted again, making the possibility to use a Multiple Objective
Genetic Algorithm.

.. toctree::
	Validation
	ObjectiveFunctions

.. note::
	**Installation instructions:**
	Instructions to install KaSim, BioNetGen, NFsim, and PISKaS are
	available in their source code webpages. Nonetheless, here you will find
	basic information to clone using git or download the software and install
	it.

	To install SLURM, you should have admin access to your infrastructure and an
	UNIX-based OS. Detailed instructions are provided here:
	:ref:`SLURM-instructions`

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
