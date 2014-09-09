Getting started guide for Windows systems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The installation of PyDSM involves two major steps: the installation of
its pre-requisites and the installation of PyDSM itself.

Note that not all version of PyDSM may be available with a Windows
installer.

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

Please, look at the PyDSM reference for further information on how to
take advantage of the code. It may also be helpful to check the
reference for information on the scientific papers that describe the
methods implemented in the package.

If you find this code useful, please consider citing such papers
in your work.

.. _pythonlibs site by Christoph Gohlke :
   http://www.lfd.uci.edu/~gohlke/pythonlibs/
.. _ActiveState Python : http://www.activestate.com/activepython
.. _Enthought Python : http://www.enthought.com/products/epd.php
.. _Python(x,y) : http://code.google.com/p/pythonxy/
.. _WinPython 2.7.x 32 bit : http://code.google.com/p/winpython/
.. _PyDSM download site : http://code.google.com/p/pydsm/downloads/list
.. _Netlib archive of prebuilt ATLAS libraries for Windows :
   http://www.netlib.org/atlas/archives/windows/
.. _ATLAS sourceforge site : http://math-atlas.sourceforge.net/
