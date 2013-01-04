Getting started guide for Linux systems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The installation of PyDSM involves two main steps. The installation of
its prerequisites and the installation of PyDSM itself.

Installation of the prerequisites
'''''''''''''''''''''''''''''''''

Installing the prerequisites with Linux is quite easy if your system is based on a modern Linux distribution like Ubuntu, Fedora, etc.

Steps that can be practiced using the distribution package manager
``````````````````````````````````````````````````````````````````

First of all, there are some steps that can usually be practiced using
the distribution package manager.

#. Make sure that a *Python 2.7* environment is installed on your
   system using the distribution package manager. For instance, in
   Ubuntu or Debian, you can use a (graphical) package manager to look
   for the packages ``python2.7`` and ``python2.7-dev``.

#. Similarly, assure that *Numpy*, *Scipy* and *Matplotlib* are
   installed on your system. Since many distribution offer Numpy,
   Scipy and Matplotlib both for Python 2 and for Python 3, assure
   that the version that you have installed matches your Python 2.7
   environment.

#. Use your package manager also to install (or verify the
   installation) of *Cython*.

#. Again, use your package manager to verify that a development
   version of the *Atlas* libraries is installed on your system. Atlas
   is a provider of Lapack and Blas (and *CBlas*, the C interface to
   Blas), which are scientific libraries for linear algebra and matrix
   manipulation.

#. Use your package manager to ensure that your system has a *C
   compiler*. Most Linux system will have the ``gcc`` C compiler
   installed by default.

#. Use your package manager to verify that your system has the
   *CVXOPT* package installed. In Ubuntu and Debian, this is called
   ``python-cvxopt``. Most likely it will not be installed by default
   and you will need to ask the package manager for its installation.

#. (Optional) Use your package manager to install the *Spyder*
   development environment.

All the above is expected to take you just a few minutes.

What to do if your distribution does not provide all of the above
.................................................................

It may be the case, that your Linux distribution does not hold all the
software above in a packaged version, or that it provides packages
with not sufficiently up to date versions. In this situation, you may
need to install some of the packages above in a less automated
fashion. Most likely, requisites about Python itself, the C compiler
and Atlas will be satisfied, and you may need to take care of
installing (newer) versions of Numpy, Scipy, Matplotlib, Cython and
CVXOPT. To this aim, you can either:

* Download the source packages  from the PyPy_ repository, unpack them
  and launch  the ``setup.py`` script with command  ``install``. It is
  advisable to use  the ``--user`` option with the  install command to
  make  a  personal  installation  that does  not  require  super-user
  permissions  and  does not  interfere  with  the Linux  distribution
  packaging system.

* Use ``easy_install`` or ``pip`` to automatize the steps above. These
  are scripts that automatically download packages from Pypi and
  install them. Again, remember to install with the ``--user`` option
  if you want to avoid any kind of interference with the distribution
  package manager.

Some notes on ATLAS
...................

The Atlas libraries provided by your distribution may not be fully
optimized for your actual hardware platform (i.e., specific CPU,
etc.). To get maximum performance, you may want to re-compile it. Your
distribution will certainly provide instructions for this. For
instance, Debian and Ubuntu tell you how to optimize the Atlas library
in the file ``/usr/share/doc/libatlas3-base/README.Debian``.

Steps that must be practiced manually
`````````````````````````````````````

One of the prerequisites of this code, namely *CVXPY* is too young to
be available in the Linux distributions or on Pypi. Thus it needs to
be installed manually. This is by no means hard, though.

This code is available at the CVXPY_ site. Please, look there for
documentation and for extended installation instructions.
Unfortunately, so far the `CVXPY download site`_ is still empty,
meaning that the code can only be downloaded with a specific tool for
accessing the `CVXPY download repository`_. Since some may find this
impractical, we are providing a development snapshot for convenience
at the `PyDSM download site`_.

Once, you have the source code of CVXPY, unzip it, enter its directory
and launch the `setup.py` file as::

  python setup.py install --user

(on some systems you may need to use something like ``python2.7``
instead of ``python``).

Installation of PyDSM itself
''''''''''''''''''''''''''''

After all the pre-requisites are satisfied, you may eventually proceed
to installing PyDSM itself. Download its source from the `PyDSM
download site`_, then unzip it and enter the source directory. Then
launch the ``setup.py`` file as::

   python setup.py install --user

(on some systems you may need to use something like ``python2.7``
instead of ``python``).

Alternatively, you may run some tests on the code before installing it
::

   python setup.py test

However, note that the tests are currently rather incomplete.

Using the code
''''''''''''''

To use PyDSM, open your Python interpreter (or the Spyder development
environment) and
::

  import pydsm

This command should not report any error. After issuing it, the PyDSM
functions should be available under the ``pydsm`` namespace.

Please, refer to the PyDSM reference for further information on how to
take advantage of the code. It may also be helpful to check this
documentation for information on how to reproduce the results in the
paper

  Sergio Callegari, Federico Bizzarri *"Output Filter Aware
  Optimization of the Noise Shaping Properties of ΔΣ Modulators via
  Semi-Definite Programming"*, IEEE Transactions on Circuits and
  systems - Part I: Regular Papers.

.. _PyPy : http://pypi.python.org/pypi
.. _CVXPY : http://www.stanford.edu/~ttinoco/cvxpy/
.. _CVXPY download site : http://code.google.com/p/cvxpy/downloads/list
.. _CVXPY download repository : http://code.google.com/p/cvxpy/source/checkout
.. _PyDSM download site : http://code.google.com/p/pydsm/downloads/list
