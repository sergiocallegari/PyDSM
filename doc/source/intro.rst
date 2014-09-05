Introduction
------------

PyDSM is a Python Delta Sigma Modulator toolbox. It contains tools for
experimenting with ΔΣ modulators. Currently, it is still relatively
small and mostly focused on the experimentation of different
techniques for the design of the modulator Noise Transfer Function
(NTF). Furthermore, the current codebase contains means to simulate a
generic modulator.

PyDSM is under development and shall be enriched with further
functionalities in a near future.

Highlights of the current version consist in:

* Some routines ported from the very well known `DELSIG toolbox
  <http://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox>`__
  for Matlab by R. Schreier
* The method for the design of psychoacoustically optimal modulators
  for audio signals proposed by Dunn and Sandler in 1997.
* The NTF design algorithm presented in:

  * Sergio Callegari, Federico Bizzarri *“Output Filter Aware
    Optimization of the Noise Shaping Properties of ΔΣ Modulators via
    Semi-Definite Programming”*, IEEE Transactions on Circuits and
    systems - Part I: Regular Papers, Vol. 60, N. 9,
    pp. 2352-2365. Sept. 2013. DOI: `10.1109/TCSI.2013.2239091
    <http://dx.doi.org/10.1109/TCSI.2013.2239091>`_. Pre-print
    available on `ArXiv <http://arxiv.org/abs/1302.3020>`__.

  * Sergio Callegari, Federico Bizzarri *“Noise Weighting in the
    Design of ΔΣ Modulators (with a Psychoacoustic Coder as an
    Example),”* IEEE Transactions on Circuits and Systems - Part II:
    Express Briefs, Vol. 60, N. 11, pp. 756-760. Nov. 2013. DOI:
    `10.1109/TCSII.2013.2281892
    <http://dx.doi.org/10.1109/TCSII.2013.2281892>`_. Pre-print
    available on `ArXiv <http://arxiv.org/abs/1309.6151>`__.

  The latter includes an optimal strategy for the design of
  psychoacoustically optimal modulators for audio signals,
  configurable to deal with different noise weightings
  (e.g. A-Weighting, F-Weigting, user supplied weightings, etc.)

If you find this code useful, *please consider citing the above papers
in your work.*
