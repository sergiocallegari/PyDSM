Dependecies on Python packages
``````````````````````````````

These are dependencies on Python packages that are available on the
`Python Package Index`_ (PyPI_).

Numpy_
    A powerful library for vector and matrix manipulation in Python.

Scipy_
    A collection of tools for scientific and engineering computing in
    Python.

Matplotlib_
    A widely used plotting library for Python.

CVXPY_
    A powerful Python package for modeling convex optimization
    problems in a natural mathematical syntax.

PICOS_
    Another is a user-friendly Python package providing a high-level
    interface for expressing optimization models, handling both convex
    and mixed-integer problems, and interfacing with various solvers.

CVXOPT_
    A Python package for convex optimization, including linear and
    quadratic programming, second-order cone programming, and related
    problems.

SCS_
    Another fast numerical solver for convex cone problems that can
    handle large-scale linear, second-order cone, and semidefinite
    programs

In addition, the following is required when *building* PyDSM from
source (currently necessary, since no binary wheels are provided):

Cython_
    A programming language designed to simplify writing C extensions
    for Python.

Using the *reference workflow*, dependency management is straightforward:
most packages are distributed in binary form as precompiled,
ready-to-install archives (*wheels*) for all major platforms. When you
install PyDSM with ``pip`` or equivalent tools, these dependencies are
automatically retrieved and installed without the need to compile them.

In the rare case where a wheel is not available for your platform,
``pip`` (or the corresponding tool) will fall back to building the
package from its source distribution (*sdist*). The
:doc:`external dependencies <external-dependencies>` required by PyDSM
are generally sufficient to build these source packages as well.

.. note::

   If you use the *conda workflow*, first create an environment and
   install the required dependencies directly with Conda. Recent
   versions of these packages are typically fine.

   Only then install PyDSM using ``pip``. This way, ``pip`` will
   detect that the dependencies are already present and will not try to
   reinstall or overwrite them. Allowing ``pip`` to manage these
   dependencies inside a Conda environment may otherwise cause
   conflicts or inconsistencies.


.. include:: _links.rst
