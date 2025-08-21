Getting Started
---------------

PyDSM is designed to run on all major platforms (Linux, Windows,
macOS, etc.). It is free software, and all its runtime prerequisites
are also free, ensuring that anyone can use it. The toolbox is
routinely built and tested by its developers on Linux and Windows, and
is tested on macOS on a more occasional basis.

The code is written in Python_ and includes some C (specifically
Cython) extensions for efficiency. A Python 3 interpreter is required
(version 3.10 or above).


The two workflows for Python package management
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before discussing the installation of PyDSM, it is worth recalling
that there are two major *workflows* for managing Python packages and
environments:

- the *reference workflow*, recommended by the Python Packaging
  Authority (PyPA_) and defined through PEPs, the formal documents by
  which all changes to the Python language and ecosystem are
  specified.

  This workflow relies on a central package index (PyPI_) and the
  standard Python installer (pip_) or compatible tools (for example,
  the package and project manager `uv`_). It is often referred to as
  “standard Python,” “PyPI-based,” or “pip-based” installation;

- the *conda workflow*, named after conda_, the package manager used
  in the Anaconda_ distribution. Conda is an open-source,
  cross-platform, and language-agnostic package manager and
  environment management system. It is particularly effective at
  distributing *binary software* and *platform-specific code*
  alongside pure Python packages. You use conda whenever you install
  Anaconda_, its lighter variant, Miniconda_, or Miniforge_.

  Anaconda is available on all three major platforms (Windows, Linux,
  and macOS). While it is not the only Python distribution, it is by
  far the most widely used. It is especially popular on platforms that
  lack a native package manager for installing Python and related
  tools — most notably Windows.

  For a package to be installed with the conda workflow, it must be
  built in a format different from the one used on PyPI_ and published
  through Anaconda_ or a compatible channel (such as conda-forge).

PyDSM follows the recommendations of the `Python Packaging User
Guide`_, the authoritative resource on packaging, publishing, and
installing Python projects using the reference workflow (although in
some areas it is still catching up). Accordingly, PyDSM is available
on PyPI_ and can be installed via pip_.

At present, PyDSM is not packaged for conda. Therefore, even if you
use Anaconda_, you must still rely on ``pip`` or a compatible tool to
install it.

.. note::

   While it is technically possible to use ``pip`` together with the
   conda workflow, there are important caveats to consider. Mixing
   packages installed with ``pip`` and ``conda`` in the same
   environment can cause version conflicts or broken dependencies,
   since the two tools resolve and manage packages differently.


PyDSM dependencies and prerequisites
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Apart from a standard Python environment, PyDSM requires a number of
prerequisites. Some are needed at *runtime*, while others are required
only at *build time* to install the package (since it is currently
distributed in source form only).

Most prerequisites are regular Python packages, which are handled
automatically during installation. In addition, there are some
*external dependencies* (non-Python or system libraries) that must be
installed manually in advance, often using platform-specific
procedures.


.. toctree::
   :maxdepth: 1

   external-dependencies
   python-dependencies


PyDSM installation
~~~~~~~~~~~~~~~~~~

Although the list of prerequisites may look extensive, managing them is
usually straightforward. Once the required external dependencies are in
place, PyDSM can be installed directly from its sources on PyPI with a
simple command:

.. code-block:: console

   pip install pydsm

.. note::

   Installing with ``pip`` is required even if you are following the
   *conda workflow*, since PyDSM is not currently available in Anaconda_.

   As mentioned in the section on :doc:`Python dependencies
   <python-dependencies>`, you should run ``pip install pydsm`` only
   *after* installing the other Python dependencies with ``conda``.
   This minimizes the risk of conflicts between ``pip`` and ``conda``.

As always with Python, working inside a *virtual environment* is
strongly recommended.


Testing the PyDSM installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To verify that PyDSM was installed correctly, start a Python REPL
session and type:

.. code-block:: pycon

   >>> import pydsm

The import should complete without errors.

PyDSM also includes a (still limited) test suite. Once the package has
been imported, the test suite can be run with:

.. code-block:: pycon

   >>> pydsm.test()

The tests should result in no errors, even if some warnings are
expected.


Additional tools to experiment with PyDSM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In addition to PyDSM itself, we recommend trying either a `Jupyter
notebook`_ or the `Spyder IDE`_ as environments for experimenting with
it. Jupyter provides an interactive workspace that combines code,
output, and rich text (including headings, equations, and more). Spyder
is a Python development environment tailored to scientific computing,
particularly convenient for users transitioning from MATLAB.

Both environments offer integrated online help, and PyDSM has been
designed so that its objects, functions, and modules include internal
documentation that integrates seamlessly with these systems.


.. include:: _links.rst
