# pleione

Pleione is a python3 package that implement methods that are common to
traditional modeling frameworks, and apply them to analyze Rule-Based Models.

Here you'll find the necessary documentation to install and use the methods in
Pleione. At the moment, Pleione parameterizes Rule-Based Models written either
in *BioNetGen* (https://www.csb.pitt.edu/Faculty/Faeder/?page_id=409) or *kappa*
language (https://www.kappalanguage.org/). Models are simulated with BNG2
(https://github.com/RuleWorld/bionetgen, PMID 27402907), NFsim
(https://github.com/RuleWorld/nfsim, PMID 26556387), KaSim
(https://github.com/Kappa-Dev/KaSim, PMID 29950016), or PISKaS
(https://github.com/DLab/PISKaS, PMID 29175206). Please contact us or write an
issue to include your favorite stochastic simulator to Pleione
(https://github.com/glucksfall/pleione/issues).

Pleione implements a Genetic Algorithm with elitism, on the contrary to
BioNetFit (https://github.com/RuleWorld/BioNetFit, PMID 26556387) that
implements a parents selection within a distribution probability that is inverse
to the rank. Nonetheless, Pleione's methods to parameterize Rule-Based Models
include both, a uniform or inverse to the rank probability to select models from
within an elite or all models.

Examples to run Pleione are located in https://github.com/glucksfall/pleione/tree/master/example
and in the python distribution wheel. The table with the U-test critical values is
located at the same folder and in subfolders.

The plan to add methods into Pleiades (https://github.com/glucksfall/pleiades)
includes a sensitivity analysis and a parameterization employing a Particle
Swarm Optimization protocol. You could write us if you wish to add methods into
pleione or aid in the development of them.
