External dependencies
`````````````````````

As previously noted, some parts of PyDSM (and notably the ΔΣ modulator
simulator) are implemented in languages other than Python for efficiency.
These components require the package to be *built*, i.e. compiled into
binary form. Since PyDSM is currently distributed on PyPI_ only as a
source package (no precompiled *wheels* are provided), building is always
required during installation.

As a result, PyDSM depends on a C compiler toolchain. The appropriate
toolchain varies by platform (Windows, Linux, or macOS), and detailed
instructions for each are provided below. Note that this is strictly a
*build-time dependency*: the compiler toolchain is required to *install*
or *upgrade* PyDSM, but not to *use* it.

.. note::

   Using a distribution such as Anaconda_ does not remove the need to
   manually install a C compiler toolchain—particularly on Windows,
   where Anaconda_ cannot provide the recommended compiler for building
   Python extensions.

.. note::

   On Linux only, there is an additional external dependency required
   both at build time and at runtime: a BLAS library together with its
   C interface wrapper, CBLAS.

   Please refer to the platform-specific instructions for the rationale
   behind this requirement and guidance on how to install it.


Linux specific instructions for external dependencies
.....................................................

To install PyDSM on Linux, the following items must be available
beforehand:

#. A Python interpreter

   - The recommended way is to use your distribution package manager.
   - Alternatively, you can install `uv`_ (via its official installer)
     and use it to manage a Python version. Another option is to use
     `pyenv`_.

#. A C compiler

   - Both GCC_ and Clang_ are supported. They can be installed with your
     distribution package manager.

#. A BLAS library with its CBLAS interface

   - On Linux (and only on Linux), two versions of the ΔΣ modulator
     simulator are provided:

     * one linked against the BLAS libraries bundled with SciPy,
     * one using an external BLAS implementation, typically provided
       by the system.

   - BLAS (Basic Linear Algebra Subprograms) is a standard API for
     low-level vector and matrix operations, while CBLAS is its C
     interface.
   - Because BLAS defines an interface rather than a specific
     implementation, multiple libraries exist, optimized for different
     hardware, performance targets, and licensing models.
   - Maintaining both versions of the simulator allows performance
     comparisons between SciPy’s internal BLAS and other system-level
     BLAS implementations (this may change in the future).
   - A suitable BLAS library can usually be installed from your
     distribution package manager.
   - Common options include:

     * `Netlib BLAS`_ — highly portable but relatively slow
     * OpenBLAS_ — widely used, often the default in Linux
       distributions and in SciPy
     * the BLAS component of the `Intel Math Kernel Library`_ —
       proprietary, usually free of charge, highly optimized (but not
       always available through the distribution package manager)



Windows specific instructions for external dependencies
.......................................................

To install PyDSM on Windows, the following items must be available
beforehand:

#. A Python interpreter

   - Windows does not include a usable Python interpreter out of the
     box, so Python must be installed manually.
   - Although Windows provides a Python launcher in the Microsoft Store,
     the Store build is often less predictable and may cause issues with
     build tools, compilers, or virtual environments. Alternative
     options are generally recommended:

     * Install Python using the official installers from python.org
     * Install `uv`_ (via its official installer) and use it to manage
       Python versions and environments
     * Use a Python distribution such as Anaconda_, or its smaller
       companion Miniconda_, or the community-driven Miniforge_

#. A C compiler

   - The recommended compiler is Microsoft Visual C++ (MSVC), available
     free of charge as part of the `Visual Studio Community Edition`_.
   - MSVC is the compiler used to build the official Python releases and
     precompiled wheels.
   - Recent Python versions require MSVC 14.3 (v143), included with
     Visual Studio 2022.
   - To some extent, MinGW_ can also work, though it is less commonly
     used for building Python extensions on Windows.


macOS specific instructions for external dependencies
.....................................................

To install PyDSM on macOS, the following items must be available
beforehand:

#. A Python interpreter

   - Recent versions of macOS do not ship with a usable Python
     interpreter, so Python must be installed manually. Common options
     include:

     * Installing Python using the official installers from python.org
     * Installing the homebrew_ package manager and then using it to
       install Python (e.g. ``brew install python``)
     * Installing `uv`_ (via its official installer) and using it to
       manage Python versions and environments
     * Using a Python distribution such as Anaconda_, its smaller
       companion Miniconda_, or the community-driven Miniforge_

#. A C compiler

   - The recommended compiler is Clang_, included with Apple’s Xcode_,
     the official integrated development environment available from the
     App Store.
   - Only the Xcode Command Line Tools are required to provide the
     compiler toolchain; the full Xcode IDE is not necessary.
   - The Command Line Tools can be installed directly from a terminal
     with:

     .. code-block:: console

        xcode-select --install



.. include:: _links.rst
