Changelog
---------

0.6.0
   First released version

0.6.1
   - Add project logo to the project source
   - Ship html documentation separately from main code
   - Add sample code to replicate the results in an ICECS 2012 paper

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

0.7.1
   - Fix computation of impulse response of filters that are already in
     FIR form.
   - Avoid direct access to numpy array data in Cython code. This is in
     preparation for future releases of numpy where direct access to
     array data is already deprecated.
   - Implement the synthesizeChebyshevNTF NTF design strategy from DELSIG.

0.7.2
   - Make codebase compatible with scipy 0.12.0
   - Make delsig module contain its reference delsig version
   - Minor fixes to the documentation

0.7.3
   - Apply fixes introduced in DELSIG 7.4
   - Make port of DELSIG functions more consistent with DELSIG
   - Provide better documentation to some functions
   - Bug fixes

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
