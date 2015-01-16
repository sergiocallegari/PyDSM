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
* The NTF design algorithm presented in [Cal13a]_ and [Cal13b]_.  The
  latter includes an optimal strategy for the design of
  psychoacoustically optimal modulators for audio signals,
  configurable to deal with different noise weightings
  (e.g. A-Weighting, F-Weigting, user supplied weightings, etc.)  If
  you find this code useful, *please consider citing the two papers
  in your work.*

PyDSM is free software and is licensed as detailed in the
:ref:`license` Section of this manual.

PyDSM is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  This is also detailed together with
the licensing information.
