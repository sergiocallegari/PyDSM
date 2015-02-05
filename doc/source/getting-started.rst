Getting Started
---------------

PyDSM is meant to run on all major platforms (Linux, Windows, Mac,
etc.). It is free software and it is written so that all its
prerequisites are also free. This means that anyone can try and use
it. PyDSM is routinely built and tested by its developers on Linux and
Windows 7. On a more occasional basis, it is tested on MacOs too.  The
code is written in Python and includes some C extensions for
efficiency reasons. Currently, it requires a `Python 2.7`_
environment. A port to Python 3 will likely happen in the near future.

Apart from a generic Python environment, PyDSM has the following
prerequisites:

Numpy_
    A powerful library that adds vector and matrix manipulation
    routines to Python. Currently tested with version 1.8.0 and
    current. May work with previous releases.

Scipy_
    A package of tools for science and engineering for
    Python. Currently tested with version 0.13.2 and current. May work
    with previous releases.

Matplotlib_
    A python plotting library. Currently tested with version
    1.3.1 and current. May work with previous releases.

CVXOPT_
    A free software package for convex optimization based on the
    Python programming language. Currently tested with version 1.1.6
    and current. May work with previous releases.

Furthermore, the following pre-requisites may exist in case one wants to build
from source (which is expected in Linux and MacOS and optional in Windows):

Cython_
    A language to write C extensions for the Python
    language. This is actually necessary only for compiling the
    code. Tested with versions 0.19.2 and current. May work with
    previous releases.

A C compiler
    This is only for compiling the
    extension modules and after that is not used anymore. In Linux,
    the C compiler can typically be installed from the distribution
    package manager. In MacOs, Xcode can be obtained from the app
    store. In Windows, just in case you want to build from source, the
    recommended compiler is Visual C++ Express 2008 (namely, the
    default compiler used for Python 2.7 itself), that is provided at
    no cost from Microsoft. Note that the 64 bit compiler needs to be
    purchased separately with the Microsoft Windows SDK for Windows 7
    and .NET Framework 3.5 SP1, also available at no cost from
    Microsoft.

In Linux only, two different versions of the simulator are built,
using two different sets of libraries for matrix manipulation (this is
to have a term of comparison for the simulator speed, but is likely to
change). Consequently, the Linux build also requires:

CBlas
    Blas is a library of routines for performing basic
    vector and matrix operations. CBlas is its edition tailored for C
    code. For its building, PyDSM requires a *development version* of
    CBlas, namely a version including header files. For running PyDSM
    on Linux, a mere CBlas *runtime* library is needed.  CBlas may be
    available from many sources including Netlib_, Atlas_, Openblas_,
    etc. On Linux use your distribution package manager to select one.
    The use of a CBLAS edition fully optimized for your CPU and system
    architecture is highly recommended.

Although the prerequisites appear to be numerous, their management is
actually quite easy.

A detailed getting started guide is provided for :doc:`Linux
<linux-getting-started>`, :doc:`Windows <windows-getting-started>` and
:doc:`MacOs <macos-getting-started>`.

In addition to the above prerequisites, we suggest to try Spyder_ as
an environment to run the code. This is a Python development
environment specifically tailored to suit the need of scientific
applications and to ease the learning path for those with an
experience in Matlab. Conveniently, it offers online help for
coding. The functions in PyDSM are internally documented to work with
this online help system.

Recent version of PyDSM may also benefit from alternative modeling
languages for convex optimization. The interested user may want to
install also the CVXPY_ and/or the PICOS_ Python packages.

.. toctree::
   :maxdepth: 2

   linux-getting-started
   windows-getting-started
   macos-getting-started

.. _Python 2.7: http://www.python.org/download/releases/2.7/
.. _Cython: http://www.cython.org/
.. _Netlib: http://www.netlib.org
.. _Atlas: http://math-atlas.sourceforge.net/
.. _Openblas: http://www.openblas.net/
.. _Numpy: http://sourceforge.net/projects/numpy/
.. _Scipy: http://sourceforge.net/projects/scipy/
.. _Matplotlib: http://matplotlib.org/
.. _CVXOPT: http://abel.ee.ucla.edu/cvxopt/
.. _Spyder: http://code.google.com/p/spyderlib/
.. _CVXPY: http://www.cvxpy.org
.. _PICOS: http://picos.zib.de/
