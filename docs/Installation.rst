Installation
============

There are two different ways to install pleione:

1. **Install pleione natively (Recommended).**

   *OR*

2. **Clone the Github repository.** If you are familiar with git, pleione can
   be cloned and the respective folder added to the python path. Further details
   are below.

.. note::
	**Need Help?**
	If you run into any problems with installation, please visit our chat room:
	https://gitter.im/glucksfall/pleiades

Option 1: Install pleione natively on your computer
---------------------------------------------------

The recommended approach is to use system tools, or install them if
necessary. To install python packages, you could use pip, or download
the package from `python package index <https://pypi.org/project/pleione/>`_.

1. **Install with system tools**

   With pip, you simple need to execute and pleione will be installed on
   ``$HOME/.local/lib/python3.6/site-packages`` folder or similar.

   .. code-block:: bash

	pip3 install pleione --user

   If you have system rights, you could install pleione for all users with

   .. code-block:: bash

	sudo -H pip3 install pleione

2. **Download from python package index**

   Alternatively, you could download the package (useful when pip fails to download
   the package because of lack of SSL libraries) and then install with pip. For instance:

   .. code-block:: bash

	wget https://files.pythonhosted.org/packages/a9/6b/99ae42d1bf3811ff909c8c0fde6a174fbef4d76fe25f0aa0edd13420d4d6/pleione-1.5.1-py3-none-any.whl
	pip3 install pleione-1.5.1-py3-none-any.whl --user

   .. note::
	**Why Python3?**:
	Pleione is intended to be used with >=python3.4 because python2.7 won't receive
	further development past 2020, including security updates.

   .. note::
	**pip, Python and Anaconda**:
	Be aware which pip you invoque. You could install pip3 with
	``sudo apt-get install python3-pip`` if you have system rights, or
	install python3 from source, and adding ``<python3 path>/bin/pip3`` to the
	path, or linking it in a directory like ``$HOME/bin`` which is commonly
	added to the path at login. Also be aware that, if you installed
	Anaconda, pip could be linked to the Anaconda specific version of pip, which
	will install pleione into Anaconda's installation folder.
	Type ``which pip`` or ``which pip3`` to find out the source of pip, and type
	``python -m site`` or ``python3 -m site`` to find out where is more likely
	pleione will be installed.

Option 2: Clone the Github repository
-------------------------------------

1. **Clone with git**

   The source code is uploaded and maintained through Github at
   `<https://github.com/networkbiolab/pleione>`_. Therefore, you could clone the
   repository locally, and then add the folder to the ``PYTHONPATH``. Beware
   that you should install the *pandas* package (`pandas`_) by any means.

   .. code-block:: bash

    git clone https://github.com/networkbiolab/pleione /opt
    echo export PYTHONPATH="\$PYTHONPATH:/opt/pleione" >> $HOME/.profile

   .. note::
	Adding the path to ``$HOME/.profile`` allows python to find the package
	installation folder after each user login. Similarly, adding the path to
	``$HOME/.bashrc`` allows python to find the package after each terminal
	invocation. Other options include setting the ``PYTHONPATH`` environmental
	variable in a sh file (see the example folder) or invoke ``python3 setup.py clean build install``
	to install Pleione as it was downloaded from the PyPI server.

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
