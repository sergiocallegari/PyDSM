Introduction
------------

PyDSM is a Python toolbox for Delta Sigma Modulators. It provides tools
for experimenting with ΔΣ modulators. At present, the library is still
relatively small and mainly focused on exploring different techniques
for designing the Noise Transfer Function (NTF). In addition, the code
includes functionality to simulate a generic modulator.

PyDSM is actively developed and will be extended with new features in
the near future.

Highlights of the current version include:

* Several routines ported from the well-known `DELSIG toolbox`_ for
  Matlab by R. Schreier.
* The method proposed by Dunn and Sandler (1997) [Dun97]_ for designing
  psychoacoustically optimal modulators for audio signals.
* The NTF design algorithms presented in [Cal13a]_ and [Cal13b]_.  The
  latter introduces an optimal strategy for designing psychoacoustically
  optimal modulators for audio signals, configurable for different noise
  weightings (e.g. A-Weighting, F-Weighting, user-supplied weightings,
  etc.).  If you find this code useful, *please consider citing these
  two papers in your work.*
* The NTF design algorithm proposed in [Nag12]_.
* The NTF design algorithm presented in [Cal15]_.

PyDSM is free software and is licensed as described in the
:ref:`license` section of this manual.

PyDSM is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranties of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  Further details are provided along
with the licensing information.

.. include:: _links.rst
