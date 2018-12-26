# pleione

Pleione is a python3 package that implement methods that are common to
traditional modeling frameworks, and apply them to analyze Rule-Based Models.

Here you'll find the necessary documentation to install and use the methods in
Pleione. At the moment, Pleione parameterizes Rule-Based Models written
either in *BioNetGen* or *kappa* language. Models are simulated with BNG2, 
NFsim, KaSim, or PISKaS. Please contact us or write an issue to include your
favorite stochastic simulator to Pleione (https://github.com/glucksfall/pleione/issues).

Pleione implements a Genetic Algorithm with elitism, on the
contrary of BioNetFit that implements a parents selection within
an distribution probability that is inverse to the rank. Nonetheless, pleione's
methods to parameterize Rule-Based Models include both, elitism and inverse methods.

The plan to add methods into Pleione includes a sensitivity analysis and a
parameterization employing a Particle Swarm Optimization protocol. You
could write us if you wish to add methods into pleione or aid in the development
of them.
