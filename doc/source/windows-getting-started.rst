Getting started guide for Windows systems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The installation of PyDSM involves two major steps: the installation of
its pre-requisites and the installation of PyDSM itself.

**PyDSM is currently only available for 32 bit Python environments in
Windows.** This does not mean that you cannot have PyDSM for a 64 bit
Python environment in windows, but merely that a pre-compiled PyDSM
for Windows is currently only built and tested with a 32 bit Python
environment. This is likely to change shortly.

Note that there is no issue in case you have a 64 bit Windows
platform, since you can easily install a 32 bit Python environment on
it.

If you really need a 64 bit PyDSM (e.g., because you want to simulate
delta sigma modulators for very long time spans in a single shot,
which may require access to a large amount of memory), please consider
using :doc:`Linux <linux-getting-started>`. Alternatively, you may
succeed in buildin all the PyDSM prerequisites and PyDSM itself with a
64 bit Python environment using a commercial compiler. However, with
this you are on your own (the `pythonlibs site by Christoph Gohlke`_
may be a good starting point). Some commercial distributions of
Python, such as `ActiveState Python`_ or `Enthought Python`_
[http://www.activestate.com/activepython ActiveState Python] or
[http://www.enthought.com/products/epd.php Enthought Python] may
support building PyDSM in a 64 bit environment on Windows
easily. However, this is also untested.

Installation of the prerequisites
'''''''''''''''''''''''''''''''''

Following the above introduction, these notes focus on 32 bit Windows
or 64 bit Windows with 32 bit Python.

Steps that can be practiced using a Python distribution
```````````````````````````````````````````````````````

First of all, note that rather than installing all the prerequisites
one by one, it can be convenient to install a Python distribution that
automatically includes most of them. We suggest using `Python(x,y)`_
or `WinPython 2.7.x 32 bit`_ that are scientific Python distribution
for 32 bit Windows XP/Vista/7 (and most likely 8 too). They include
Numpy, Scipy, Matplotlib, Cython, Spyder, the MinGW C compiler and
CVXOPT. Obviously, it is also possible to use other Python
distributions, or to install the prerequisites one by one, but with
this you are on your own.

Note that although Python(x,y) and WinPython include most of the
prerequisites of PyDSM, it may not install all of them by
default. Thus, during the setup of PyDSM make sure that Numpy, Scipy,
Matplotlib, Cython, Spyder, and CVXOPT are actually present on your
system.

Setting up a Python distribution such as Python(x,y) or WinPython is
expected to take just a few minutes.

Steps that must be practiced manually
`````````````````````````````````````

One of the prerequisites of PyDSM, namely *CVXPY*, is too young to be
available in Python distributions. Thus, it needs to be installed
manually. This is by no means hard, though.

This code is available at the CVXPY_ site. Please, look there for
documentation and for extended installation instructions.
Unfortunately, so far the `CVXPY download site`_ is still empty,
meaning that the code can only be downloaded with a specific tool for
accessing the `CVXPY download repository`_. Since some may find this
impractical, we are providing a development snapshot for convenience
at the `PyDSM download site`_. This is available either in source form
or as an installer (easier).

In case you want to use the installer, just download the CVXPY
``.exe`` file and run it. Some Python distributions such as WinPython
let you just drop the file on a *package manager* window to install
it.

Otherwise, once you have the source code of CVXPY, unzip it, enter its
directory and launch the `setup.py` file as::

  python setup.py install --user

The ``--user`` option is there to make a personal installation.

Warning
.......

There is a bug in the Python ``setuptools`` that *may* cause Spyder
(and possibly other Python programs) not to start after the use of
setup.py illustrated above. Thus, after the installation of CVXPY,
check if Spyder works. If it does not start, do the following:

#. Go to the directory where your personal Python packages get
   installed. This may be something like ``C:\Users\<your
   username>\AppData\Roaming\Python\Python2.7\site-packages``

#. Move away the ``setuptools.pth`` file (e.g. rename it to
   ``setuptools.pth.bak``) since it messes up the directives by which
   the Python environment searches its components.


Installation of PyDSM itself
''''''''''''''''''''''''''''

After all the prerequisites above are satisfied, you may eventually
proceed to installing PyDSM itself. Note that, differently from Linux,
in Windows there is no dependency on an external CBlas library.

Since version 0.7.0 an installer is available, simplifying the setup
on Windows. To use it, just download the pydsm ``.exe`` file from the
`PyDSM download site`_, and run it (or drop it on your Python
distribution *package manager* if you are using a Python distribution
providing this facility).

Alternatively, if you prefer compiling the source code, download it
from the `PyDSM download site`_, then unzip it. Finally, open a
command prompt inside the PyDSM source folder and launch the
``setup.py`` file as::

   python setup.py install --user

One of the advantages of using the source version is that you may run
some tests on the code before installing it ::

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
take advantage of the code. It may also be helpful to check the
reference for information on how to reproduce the results in the papers

  Sergio Callegari, Federico Bizzarri *"Output Filter Aware
  Optimization of the Noise Shaping Properties of ΔΣ Modulators via
  Semi-Definite Programming,"* IEEE Transactions on Circuits and
  systems - Part I: Regular Papers.

  Sergio Callegari, Federico Bizzarri *"Should ΔΣ modulators used in
  AC motor drives be adapted to the mechanical load of the motor?,"*
  Proceedings of the 19th IEEE International Conference on
  Electronics, Circuits and Systems (ICECS), 2012, pp. 849 - 852.

  Sergio Callegari, Federico Bizzarri *"Noise Weighting in the
  Design of ΔΣ Modulators (with a Psychoacoustic Coder as an
  Example),"* IEEE Transactions on Circuits and Systems - Part II:
  Express Briefs. To appear in 2013.

If you find this code useful, please consider citing the above papers
in your work.

.. _pythonlibs site by Christoph Gohlke :
   http://www.lfd.uci.edu/~gohlke/pythonlibs/
.. _ActiveState Python : http://www.activestate.com/activepython
.. _Enthought Python : http://www.enthought.com/products/epd.php
.. _Python(x,y) : http://code.google.com/p/pythonxy/
.. _WinPython 2.7.x 32 bit : http://code.google.com/p/winpython/
.. _CVXPY : http://www.stanford.edu/~ttinoco/cvxpy/
.. _CVXPY download site : http://code.google.com/p/cvxpy/downloads/list
.. _CVXPY download repository : http://code.google.com/p/cvxpy/source/checkout
.. _PyDSM download site : http://code.google.com/p/pydsm/downloads/list
.. _Netlib archive of prebuilt ATLAS libraries for Windows :
   http://www.netlib.org/atlas/archives/windows/
.. _ATLAS sourceforge site : http://math-atlas.sourceforge.net/
