Getting started guide for Windows systems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The installation of PyDSM involves two major steps: the installation
of its :doc:`pre-requisites <getting-started>` and the installation of
PyDSM itself.

Note that *a windows installer may not be provided for all version of
PyDSM*.

Installation of the prerequisites
'''''''''''''''''''''''''''''''''

Steps that can be practiced using a Python distribution
```````````````````````````````````````````````````````

First of all, note that rather than installing all the prerequisites
one by one, it can be convenient to install a Python distribution that
automatically includes most of them. We suggest using `Python(x,y)`_
or `WinPython`_ that are scientific Python distribution for Windows
Vista/7/8. WinPython is also available in a flavor for 64 bit Windows.
They include Numpy, Scipy, Matplotlib, Cython, Spyder. Obviously, it
is also possible to use other Python distributions, or to install the
prerequisites one by one, but with this you are on your own.

Note that although Python(x,y) and WinPython include most of the
prerequisites of PyDSM, they may not install all of them by
default. Thus, during the setup of PyDSM make sure that all the PyDSM
prerequisites that are available in the python distribution are also
flagged for installation.

Setting up a Python distribution such as Python(x,y) or WinPython is
expected to take just a few minutes.

Steps that may not be possible with a Python distribution
`````````````````````````````````````````````````````````

It may be the case that the adopted python distribution misses some
required packages. Actually, most of the packages required by PyDSM
are quite common such as `Numpy`_, `Scipy`_, and `Matplotlib`_ and should
definitely be available in any scientific oriented
distribution. The only less common pre-requisite is `CVXOPT`_.

In any case, if something is missing, it needs to be manually
installed. Typically, this will not require any compilation, since
there are distributions of pre-built python packages for Windows. A
good starting point can be looking at (the `pythonlibs site by
Christoph Gohlke`_).  Remember to pick the installers that are correct
for your environment (32 or 64 bit). Also, try avoid mixing packages
from different sources, since when one deals with non-pure python
packages (namely, when packages that take advantage of C or C++
libraries), subtle incompatibilities may apper. For the same reason,
assure that the Python distribution that you are using and the package
installers are compatible with each other.

A word of caution
`````````````````

Currently, the state of affairs with respect to free Python
distributions for Windows is not completely ideal.  A few
distributions now ship without `CVXOPT`_ (notably `WinPython`_ has
dropped this package at version 2.7.6.3) and at the same time the
`CVXOPT`_ 1.1.7 package provided at the `pythonlibs site by
Christoph Gohlke`_ seems to have problems with them. Notably there are
issues when using `CVXOPT`_ 1.1.7 from the `pythonlibs site by
Christoph Gohlke`_ and `WinPython`_ the latest `WinPython`_ 2.7.6.4.

The suggestion is to use `WinPython`_ 2.7.6.2 (the latest version
shipping with `CVXOPT`_).


Installation of PyDSM itself
''''''''''''''''''''''''''''

After all the prerequisites above are satisfied, you may eventually
proceed to installing PyDSM itself. Note that, differently from Linux,
in Windows there is no dependency on an external CBlas library.

Since version 0.7.0 an installer is available, simplifying the setup
on Windows. Since version 0.9.1 the installer is also available for 64
bit Windows. To use it, just download the pydsm ``.exe`` file from the
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
.. _Python(x,y) : http://code.google.com/p/pythonxy/
.. _WinPython : http://code.google.com/p/winpython/
.. _PyDSM download site : https://code.google.com/p/pydsm/wiki/download?tm=2
.. _Numpy: http://sourceforge.net/projects/numpy/
.. _Scipy: http://sourceforge.net/projects/scipy/
.. _Matplotlib: http://matplotlib.org/
.. _CVXOPT: http://abel.ee.ucla.edu/cvxopt/
