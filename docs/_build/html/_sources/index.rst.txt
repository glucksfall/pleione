Welcome to pleione's documentation!
===================================

Pleione is a python3 package that implement methods that are common to
traditional modeling frameworks, and apply them to analyze Rule-Based Models.

Here you'll find the necessary documentation to install and use the methods in
pleione. At the moment, pleione parameterizes Rule-Based Models written
either in *BioNetGen* (`BioNetGen`_) or *kappa* language (`Kappa`_). Models are
simulated with BNG2 (`BioNetGen2`_), NFsim (`NFsim`_), KaSim (`KaSim`_), or
PISKaS (`PISKaS`_). Pleione implements a Genetic Algorithm with elitism, on the
contrary of BioNetFit (`BioNetFit`_) that implements a parental selection within
an inverse distribution. Nonetheless, pleione's methods to parameterize
Rule-Based Models include both, elitism and inverse methods.

The plan to add methods into pleione includes a sensitivity analysis and a
parameterization method employing a Particle Swarm Optimization protocol. You
could write us if you wish to add methods into pleione or aid in the development
of them.

.. toctree::
   :maxdepth: 3

   Installation
   ParameterEstimation
   Python3
   SLURM


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

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
