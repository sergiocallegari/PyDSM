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

* Some routines ported from the very well known DELSIG Toolbox for
  Matlab by R. Schreier
* The method for the design of psychoacoustically optimal modulators
  for audio signals proposed by Dunn and Sandler in 1997.
* The NTF design algorithm presented in:

    Sergio Callegari, Federico Bizzarri *"Output Filter Aware
    Optimization of the Noise Shaping Properties of ΔΣ Modulators via
    Semi-Definite Programming,"* IEEE Transactions on Circuits and
    Systems - Part I: Regular Papers. To appear in 2013.

    Sergio Callegari, Federico Bizzarri *"Noise Weighting in the
    Design of ΔΣ Modulators (with a Psychoacoustic Coder as an
    Example),"* IEEE Transactions on Circuits and Systems - Part II:
    Express Briefs. To appear in 2013.

  The latter includes an optimal strategy for the design of
  psychoacoustically optimal modulators for audio signals,
  configurable to deal with different noise weightings
  (e.g. A-Weighting, F-Weigting, user supplied weightings, etc.)

If you find this code useful, please consider citing the above papers
in your work.
