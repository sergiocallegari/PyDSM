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
