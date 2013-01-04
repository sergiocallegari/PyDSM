Getting Started
---------------

This software is meant to run on all major platforms (Linux, Windows,
Mac, etc.). It is free software written so that all its prerequisites
are free too. This means that anyone can try and use it.

The code is written in Python with some C extensions for efficiency
and requires a `Python 2.7`_  environment. A port to Python
3 will likely happen in the near future.

Apart from the generic Python environment, there are other
prerequisites too:

Cython_ :
    A language to write C extensions for the Python language. This is
    actually necessary only for compiling the code. Tested with
    versions 0.16.x and 0.17.x, should work with later releases too.

A development version of the CBlas library :
    Blas is a library of routines for performing basic vector and
    matrix operations. CBlas is its edition suitable for C code
    development. By development version, version of Cblas incluing
    headers files is intended. The latter are used to compile the code
    only, while the Cblas runtime is necessary all the time. CBlas may
    be available from many sources including Netlib_, etc.

A C compiler :
    This is used only for compiling the code.

Numpy_ :
    A powerful library that adds vector and matrix manipulation
    routines to Python. Currently tested with version 1.6.2, should
    work with other releases too.

Scipy_ :
    A package of tools for science and engineering for
    Python. Currently tested with version 0.11. Should work with later
    releases too.

Matplotlib_ :
    A python plotting library. Currently tested with version 1.1.1,
    should work with later releases too.

CVXOPT_ :
   A free software package for convex optimization based on the Python
   programming language. Currently tested with version 1.1.5.

CVXPY_ :
    A free software package for modeling convex optimization problems
    in Python. Currently tested with version 0.0.1 of the library,
    should work with later releases too.

Although the prerequisites appear to be numerous, their management is
actually quite easy.

A detailed getting started guide is provided for :doc:`Linux
<linux-getting-started>` and Windows.

In addition to the above prerequisites, we suggest to try Spyder_ as
an environment to run the code. This is a Python development
environment specifically tailored to suit the need of scientific
applications and to ease the learning path for those with an
experience in Matlab. It offers online help for coding. The functions
in PyDSM are internally documented to work with this online help
system.

.. toctree::
   :maxdepth: 2

   linux-getting-started
   windows-getting-started

.. _Python 2.7: http://www.python.org/download/releases/2.7/
.. _Cython: http://www.cython.org/
.. _Netlib: http://www.netlib.org
.. _Numpy: http://sourceforge.net/projects/numpy/
.. _Scipy: http://sourceforge.net/projects/scipy/
.. _Matplotlib: http://matplotlib.org/
.. _CVXOPT: http://abel.ee.ucla.edu/cvxopt/
.. _CVXPY: http://www.stanford.edu/~ttinoco/cvxpy/
.. _Spyder: http://code.google.com/p/spyderlib/
