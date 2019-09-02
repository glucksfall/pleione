Add a fitness function to Pleione
=================================

Each simulator are provided with two scripts that calculate errors. They are
located at the same path as the main scripts that calibrate. Inside each,
there is a template intended with instructions:

   .. code-block:: bash

	# Fitness Calculation Template:
	if set(args.error).issuperset(set(['the-acronysm'])):
		1. func = 0

		2. func = an algebraic expression combining the data average (data_avrg), data standard deviation (data_stdv), simulation average (sims_stdv),
		simulation standard deviation (sims_stdv), single experimental files (data.loc[i]), and/or simulation files (sims.loc[j])
		Note1: Perform two for-loops if using data.loc[i] and sims.loc[j].
		Note2: Please consider these variables are DataFrames, meaning that multiplication and division are methods (e.g. df1.division(df2))

		3. Drop NaN values (from experimental time points without simulated values, or simulated values without experimental data)
		with dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all'). Also transform Inf values with replace([numpy.inf, -numpy.inf], numpy.nan)

		4. Sum the two dimensions, and return a 6 float points scientific notation number (0 float points for statistical tests):
		error['the-acronysm'] = '{:.6e}'.format(func.dropna(axis = 0, how = 'all').dropna(axis = 1, how = 'all').sum().sum())

To use:

1) Define an acronysm for your fitness function and replace ``the-acronysm`` in the template

2) Define func as an operation of DataFrames: data_avrg, data_stdv, sims_avrg, sims_stdv, data.loc[i], and sims.loc[j]

.. note::
	Do not use data.iloc[i] or sims.iloc[i] as they provide wrong access to the data structure, possibly returning wrong calculations.
