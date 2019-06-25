Welcome to Pleione's documentation!
===================================

Pleione is a python3 package that implement methods that are common to
traditional modeling frameworks, and apply them to analyze Rule-Based Models.

Here you'll find the necessary documentation to install and use the methods in
Pleione. At the moment, Pleione parameterizes Rule-Based Models written
either in *BioNetGen* (`BioNetGen`_) or *kappa* language (`Kappa`_). Models are
simulated with BNG2 (`BioNetGen2`_, PMID `27402907`_), NFsim (`NFsim`_, PMID `26556387`_), KaSim (`KaSim`_, PMID `29950016`_), or
PISKaS (`PISKaS`_, PMID `29175206`_). Please contact us or write an issue to include your
favorite stochastic simulator to Pleione (https://github.com/glucksfall/pleione/issues).

Pleione implements a Genetic Algorithm with elitism, on the
contrary to BioNetFit (`BioNetFit`_, PMID `26556387`_) that implements a parents selection within
a distribution probability that is inverse to the rank. Nonetheless, Pleione's
methods to parameterize Rule-Based Models include both, a uniform or inverse to the rank
probability to select models from within an elite or all models.

The plan to add methods into Pleiades (`pleiades`_) includes a sensitivity analysis and a
parameterization employing a Particle Swarm Optimization protocol. You
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

.. _27402907: https://www.ncbi.nlm.nih.gov/pubmed/27402907
.. _26556387: https://www.ncbi.nlm.nih.gov/pubmed/26556387
.. _29950016: https://www.ncbi.nlm.nih.gov/pubmed/29950016
.. _29175206: https://www.ncbi.nlm.nih.gov/pubmed/29175206
.. _26556387: https://www.ncbi.nlm.nih.gov/pubmed/26556387

.. _pleiades: https://github.com/glucksfall/pleiades
