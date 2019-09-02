Algebraic Objective Functions
=============================

Here are the formulas to calculate the error between multiple data sets and
simulations files.

.. hlist::
	:columns: 2

	* | **Squared Difference of two Averages (SDA; formerly Mean Square Error, MSE):**
	  |
	  |  :math:`\left( \frac{1}{m} \sum_{i=1}^{m} \mathrm{exp}_i - \frac{1}{n} \sum_{j=1}^{n} \mathrm{sim}_j \right) ^{2}`

	* | **Absolute Difference of two Averages (ADA; formerly Mean Absolute Error, MAE):**
	  |
	  |  :math:`\mathrm{abs} \left( \frac{1}{m} \sum_{i=1}^{m} \mathrm{exp}_i - \frac{1}{n} \sum_{j=1}^{n} \mathrm{sim}_j \right)`

	* | **Pair-Wise Square Deviation (PWSD):**
	  |
	  |  :math:`\frac{1}{mn} \sum_{i=1}^{m} \sum_{j=1}^{n} \left({\mathrm{exp}_i - \mathrm{sim}_j } \right)^{2}`

	* | **Absolute Pair-Wise Deviation (APWSD):**
	  |
	  |  :math:`\frac{1}{mn} \sum_{i=1}^{m} \sum_{j=1}^{n} \mathrm{abs} \left( \mathrm{exp}_i - \mathrm{sim}_j \right)`

	* | **Normalized Pair-Wise Square Deviation (NPWSD):**
	  |
	  |  :math:`\frac{1}{mn} \sum_{i=1}^{m} \sum_{j=1}^{n} \left( \frac{ \mathrm{exp}_i - \mathrm{sim}_j }{ \mathrm{exp}_i } \right)^{2}`

	* | **Absolute Normalized Pair-Wise Deviation (ANPWSD):**
	  |
	  |  :math:`\frac{1}{mn} \sum_{i=1}^{m} \sum_{j=1}^{n} \mathrm{abs} \left( \frac{ \mathrm{exp}_i - \mathrm{sim}_j }{ \mathrm{exp}_i } \right)`

	* | **Sum of SQuares (SSQ):**
	  |
	  |  :math:`\sum_{i=1}^{m} \sum_{j=1}^{n} \left({\mathrm{exp}_i - \mathrm{sim}_j } \right)^{2}`

	* | **Chi-Square (CHISQ):**
	  |
	  |  :math:`\sum_{i=1}^{m} \sum_{j=1}^{n} \left( \frac{ \mathrm{exp}_i - \mathrm{sim}_j }{ \sigma_{\mathrm{exp}} } \right)^{2}`

	* | **Mean Normalized Square Error (MNSE):**
	  |
	  |  :math:`\sum_{i=1}^{m} \sum_{j=1}^{n} \left( \frac{ \mathrm{exp}_i - \mathrm{sim}_j }{ \overline{\mathrm{exp}} } \right)^{2}`

.. note::
	**Need a different Objective Function?** The code that calculates the error
	is separated from the main Genetic Algorithm. This make useful to encode
	other Objective Functions if the already implemented does not apply to your
	necessities. You could contact us to add your function to the pleione
	package. See the ``fitness.py`` file and add yours favourite function.

	Reach us on GitHub to add yours fitness function to the code, as it might help
	other users.
