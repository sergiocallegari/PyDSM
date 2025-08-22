Changelog
---------

0.15.2
   - Update the documentation and prepare it for readthedocs

0.15.1
   - Link simulator extension with ``-lblas``, rather than ``-lcblas``
     on Linux, to improve compatibility with recent distros

0.15.0
   - build process relies on modern Python practices (PEP 517);
   - support for modern Python (>=3.10);
   - extension modules now built in Python 3 mode;

0.14.0.0
   - Fix compatibility with cvxpy 1.0.x
   - Fix compatibilit with PICOS 1.2.0
   - Fix compatibility with matplotlib 3
   - Fix compatibility with sphinx 2.2.0 (html docs only)
   - Minor bugfixes

0.13.0.1
   - Fix compatibility with most recent cvxpy
   - Minor bugfixes

0.13.0.0
   - Python 3 compatibility!
   - Bugfixes

0.12.0.1
   - Fix crashes due to incorrect management of parameters in deprecated
     functions

0.12.0.0
   - Switch to new version numbering scheme also in view of
     PEP 440
   - Install as zipped package
   - Modify management of optional parameters in some functions.
     Note that this implies a minor API breakage.
   - Implement new hybrid NTF design method
   - Slightly improve accuracy of weighted NTF design functions
   - Improve quantization_noise_gain function
   - Improve html documentation
   - Let multiple modelers be selected in functions using convex
     optimization (cvxpy, cvxpy_old and picos are supported)
   - Implement some more functions in delsig module (axisLabels, rms)
   - Many small fixes

0.11.0
   - Switch to setuptools for building
   - Improve testing framework
   - Better management of optional parameters in some functions
   - API cleanups
     (a best effort has been put in retaining back-compatibility. Some
     back compatibility functions are deprecated and will be removed)
   - Code style improvements

0.10.1
   - Fix setup script for compatibility with MacOs
   - Provide getting started guide for MacOs
   - Minor fixes in documentation

0.10.0
   - Ready setup script for distribution
   - Implement ``partitionABCD`` in delsig module
   - Make delsig module PEP8 compliant
   - Some minor improvements to utility functions

0.9.1
   - Apply some fixes to the modulator simulator
   - Make building for 64 bit windows possible
   - Implement ``clans`` NTF design method
   - Implement ``minmax`` NTF design method
     (only single band LP, so far)

0.9.0
   - Include a local version of the discontinued cvxpy package
     by Tomas Tinoco de Ribera. This is a temporary measure

0.8.3
   - Fix a typo in the fast DSM simulator
     (only affecting case where modulator structure is passed in ABCD form)
   - Add example from ICECS 2013 paper
   - Enhance ``quantization_weighted_noise_gain`` function

0.8.2
   - Fix some licensing issues

0.8.1
   - Prevent synthesizeNTF failure if there are no zeros to optimize
   - Improve some docstrings
   - Remove some spurious imports
   - Fix normalization in quantization_noise_gain functions

0.8.0
   - Add NTF design method based on a noise weighting function
   - Provide a new module with standard audio weighting functions
   - Provide a new module with ISO 226 equal loudness contours
   - Provide a new module with NTF design methods for psychoacoustically
     optimal modulators for audio signals
   - Fix a regression in ``ds_optzeros`` introduced with version 0.7.3
     and preventing some example code from running
   - Add new examples from a recently published TCAS-II paper
   - Use ``'ba'`` specifier for requiring filters in numerator/denominator form
   - Make ``evalTF`` function more robust against complex overflow
   - Bug fixes

0.7.3
   - Apply fixes introduced in DELSIG 7.4
   - Make port of DELSIG functions more consistent with DELSIG
   - Provide better documentation to some functions
   - Bug fixes

0.7.2
   - Make codebase compatible with scipy 0.12.0
   - Make delsig module contain its reference delsig version
   - Minor fixes to the documentation

0.7.1
   - Fix computation of impulse response of filters that are already in
     FIR form.
   - Avoid direct access to numpy array data in Cython code. This is in
     preparation for future releases of numpy where direct access to
     array data is already deprecated.
   - Implement the synthesizeChebyshevNTF NTF design strategy from DELSIG.

0.7.0
   - Dropped dependency on ATLAS on Windows. Now using the blas functions
     made available via scipy. The linux version still uses ATLAS that has
     a little performance advantage.
   - Much simpler installation on Windows
   - API changes: renamed ``synthezize_ntf_from_filter_ir`` into
     ``synthezize_ntf_from_filter_imp``; swapped param order in
     ``q0_from_filter_imp_response``.
   - Fixed passing of options to ``synthesize_ntf_from_q0`` and
     ``synthesize_ntf_from_filter_imp``

0.6.1
   - Add project logo to the project source
   - Ship html documentation separately from main code
   - Add sample code to replicate the results in an ICECS 2012 paper

0.6.0
   First released version
