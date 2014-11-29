Getting started guide for MacOs systems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The installation of PyDSM involves two main steps. The installation of
its prerequisites :doc:`pre-requisites <getting-started>` and the
installation of PyDSM itself.


Installation of the prerequisites
'''''''''''''''''''''''''''''''''

Installing the prerequisites with MacOs is quite easy since most of
them are automatically provided by the platform.


Installation of the Python environment
``````````````````````````````````````

A Python 2.7 environment is ready available on all recent MacOS
versions (from Lion on). Still, it is strongly recommended to install
a stand alone Python distribution (such as Anaconda_) to minimize the
risk of interference with the system. The advantage of using a Python
distribution is that it also provides most of the packages that PyDSM
needs (e.g., Numpy, Scipy, Matplotlib, etc.). In case you do not want
to use a Python distribution, you will first need to install a C
compiler (see the next section) and then to install all the individual
pre-requisites. The ``pip`` command can do it for you.


Installation of the C compiler
``````````````````````````````

To make sure that your system includes a C compiler, you can install
Xcode for the Mac App Store.


Installation of PyDSM itself
''''''''''''''''''''''''''''

After the Python environment and the C compiler are set up, you may
eventually proceed to installing PyDSM itself. As of today, the
package is available from the PyPi_ repository, so that it is
sufficient to issue the following command:

   pip install pydsm

This downloads, unpacks and installs the package in a single step. In
case there are pre-requisites that need to be satisfied, ``pip`` also
tries to download and install them.  The ``--user`` option can be
conveniently added to the command above to ensure a personal
installation that does not require administrator permissions and does
not interfere with the Linux distribution packaging system.

As an alternative, you can manually download the package either from
PyPi_ (recommended) of from the `PyDSM download site`_. After the
download, you need to expand the archive and launch the
``setup.py`` file as::

   python setup.py install

As before, you may want to provide the ``--user`` option to make a
personal installation.


Testing the code
''''''''''''''''

PyDSM includes a (rather incomplete for the moment) set of self tests.
These can be run from the source tree as by issuing the command::

  python setup.py test

Alternatively, once the package is installed, the tests can be run by
opening a python interpreter and by typing::

  import pydsm
  pydsm.test()


Using the code
''''''''''''''

To use PyDSM, open your Python interpreter (or the Spyder development
environment) and
::

  import pydsm

This command should not report any error. After issuing it, the PyDSM
functions should be available under the ``pydsm`` namespace.

Please, look at the PyDSM reference for further information on how to
take advantage of the code. It may also be helpful to check the
reference for information on the scientific papers that describe the
methods implemented in the package.

If you find this code useful, please consider citing such papers
in your work.

.. _Anaconda : https://store.continuum.io/cshop/anaconda/
.. _CVXOPT: http://abel.ee.ucla.edu/cvxopt/
.. _PyPi : http://pypi.python.org/pypi
.. _PyDSM download site : https://code.google.com/p/pydsm/wiki/download?tm=2
