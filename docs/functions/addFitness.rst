Add a fitness function to Pleione
=================================

Each simulator are provided with two scripts that calculate errors. They are
located at the same path as the principal scripts. Inside, they have a template
intended with instructions:

   .. code-block:: bash

	# Fitness Calculation Template:
	if set(args.error).issuperset(set(['the-acronysm'])):
		func = 0
		func = an algebraic expression combining the data average (data_avrg), data variance (data_stdv), simulation average (sims_stdv),
		single experimental files (data.loc[i]) and/or simulation files (sims.loc[i]).
		# Please consider this variables are DataFrames, meaning that division is a method (pandas.DataFrame.division)
		# Please calculate average or standard deviation values from data.loc[i] and sims.loc[i] if they are needed from them (as in MSE)

		error['acronysm'] = '{:.6e}'.format(func.dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all').sum().sum())
		# drop NaN values (from experimental data without simulation point or vice-versa), sum the two dimensions, and return a 6 float points scientific notation number

To use:

1) Define an acronysm for your fitness function and replace "the-acronysm"

2) Define func as an operation of DataFrames: data_avrg, data_stdv, sims_stdv, data.loc[i], and sims.loc[i]

.. note::
	simulator-doerror.py scripts calculates one single fitness function at the time.
	The Mean Square Error has code to calculate the average from data and simulations.
