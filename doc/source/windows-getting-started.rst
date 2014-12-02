Getting started guide for Windows systems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The installation of PyDSM involves two major steps: the installation
of its :doc:`pre-requisites <getting-started>` and the installation of
PyDSM itself.

It is worth anticipating that building software from source in Windows
is somehow less comfortable than on the other platforms where PyDSM is
supported (this is particularly true of 64bit Windows). For this
reason *binary installers* may be available for Windows as a
convenience. Due to limited resources, not all versions of PyDSM may
ship a binary installer, though.  Binary installers make the
installation simpler, since it is not necessary to have a C
development environment. However, they may also create problems when
application binary interfaces are mismatched.


Installation of the prerequisites
'''''''''''''''''''''''''''''''''


Steps that can be practiced using a Python distribution
```````````````````````````````````````````````````````

First of all, note that rather than installing all the prerequisites
one by one, it is convenient to install a Python distribution that
automatically includes most of them. The recommendation is to use
`Python(x,y)`_ or `WinPython`_ that are scientific Python
distributions for Windows Vista/7/8. WinPython is also available in a
flavor for 64 bit Windows.  They include Numpy, Scipy, Matplotlib,
Cython, Spyder. Obviously, it is also possible to use other Python
distributions, or to install the prerequisites one by one, but with
this you are on your own.

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
are quite common such as Numpy, Scipy, and Matplotlib and should
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
bit Windows. To use it, just download the pydsm ``.exe`` file from
Pypi_ (recommended) or the `PyDSM download site`_, and run it (or drop
it on your Python distribution *package manager* if you are using a
Python distribution providing this facility).


Building from source
''''''''''''''''''''

To build from source, you need a C compiler. Unfortunately, you may
not need just a C compiler. In fact a *specific* C compiler may be
required to avoid subtle inconsistencies in binary interfaces. The
recommended compiler for working with Python 2.7 is `Microsoft Visual
C++ Express 2008`_. This is the same compiler used to build the
official Python 2.7 interpreter. You can get the compiler at no cost
from Microsoft, using the previous link. Unfortunately, the 2008
edition of Visual C++ Express does not contain any 64 bit compiler,
that you need if you are working on a 64 bit edition of Windows and
you deploy a 64 bit Python distribution. You will need to install the
64 bit compiler separately. This is available with the `Microsoft
Windows SDK for Windows 7 and .NET Framework 3.5 SP1`_.  Only very
recently, Microsoft has started providing a `Microsoft Visual C++
Compiler for Python 2.7`_ which may contain all that is needed in a
single package (but still there has not been any occasion to test it
with PyDSM - your feedback is welcome). Very reassuringly Microsoft
says "This compiler package is entirely unsupported".

Once you have a working C compiler, you can install PyDSM as simply as
issuing::

    pip install pydsm

In the above command, you may want to add the ``--user`` option to
force a personal installation that does not require administrator
privilege and has a lower risk of interfering with your platform
Python installation.

As an alternative, you can manually download the package either from
PyPi_ (recommended) of from the `PyDSM download site`_. After the
download, you need to expand the archive and launch the
``setup.py`` file as::

   python setup.py install

As before, you may want to provide the ``--user`` option to make a
personal installation.


So why cannot I build all my pre-requisites from source?
````````````````````````````````````````````````````````

If you have been patient enough to install the Visual C++ compiler,
you may wonder why not installing all the PyDSM pre-requisites from
source, including CVXOPT, in order to avoid any possible application
binary interface inconsistency.

The quick answer is *because you do not have a development version of
the BLAS libraries*. These are matrix manipulation libraries that are
at the basis of the numerical computations required by Numpy, CVXOPT
and PyDSM. Unfortunately, a distribution of the BLAS libraries for
Windows that is both free and easy to use does not exist. On one side
you have very expensive packages such as the Intel Math Kernel
Library.  On the other side you have free BLAS editions that are
either very inefficient (e.g., the reference Netlib_ implementation) or
relatively efficient, but a pain to compile on Windows
(e.g., Atlas_). Possibly, the advent of Openblas_, that provides
pre-built binaries for windows may change the scenario in the future.


Testing the code
''''''''''''''''

PyDSM includes a (rather incomplete for the moment) set of self tests.
Once the package is installed, the tests can be run by
opening a python interpreter and by typing::

  import pydsm
  pydsm.test()

Alternatively, when building from the source tree (see above), the
tests can be run by issuing the command::

  python setup.py test


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
.. _PyPi : http://pypi.python.org/pypi
.. _PyDSM download site : https://code.google.com/p/pydsm/wiki/download?tm=2
.. _CVXOPT: http://abel.ee.ucla.edu/cvxopt/
.. _Microsoft Visual C++ Express 2008 : http://go.microsoft.com/?linkid=7729279
.. _Microsoft Visual C++ Compiler for Python 2.7 : http://www.microsoft.com/en-us/download/details.aspx?id=44266
.. _Microsoft Windows SDK for Windows 7 and .NET Framework 3.5 SP1 : http://www.microsoft.com/en-us/download/details.aspx?id=3138
.. _Netlib: http://www.netlib.org
.. _Atlas: http://math-atlas.sourceforge.net/
.. _Openblas: http://www.openblas.net/
