PyDSM
=====

This is a Python/Scipy toolbox for the design and simulation of
ΔΣ Modulators, with emphasis on digital modulators.

Currently, it focuses on tools for the design of the modulator NTF. It
also includes a fast simulator for digital modulators. Furthermore, it
includes a Python/Scipy port of some functions from the `DELSIG
toolbox
<http://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox>`__
by R. Schreier.

As a highlight, the code includes an original NTF design technique
fully described in the papers:

* Sergio Callegari, Federico Bizzarri *“Output Filter Aware
  Optimization of the Noise Shaping Properties of ΔΣ Modulators via
  Semi-Definite Programming”*, IEEE Transactions on Circuits and
  systems - Part I: Regular Papers, Vol. 60, N. 9,
  pp. 2352-2365. Sept. 2013. DOI: `10.1109/TCSI.2013.2239091
  <http://dx.doi.org/10.1109/TCSI.2013.2239091>`_. Pre-print available
  on `ArXiv <http://arxiv.org/abs/1302.3020>`__.
* Sergio Callegari, Federico Bizzarri *“Noise Weighting in the Design
  of ΔΣ Modulators (with a Psychoacoustic Coder as an Example),”* IEEE
  Transactions on Circuits and Systems - Part II: Express Briefs,
  Vol. 60, N. 11, pp. 756-760. Nov. 2013. DOI:
  `10.1109/TCSII.2013.2281892
  <http://dx.doi.org/10.1109/TCSII.2013.2281892>`_. Pre-print available
  on `ArXiv <http://arxiv.org/abs/1309.6151>`__.

If you find the code useful, *please cite these papers in your work*.

----

The code is available on the Python Package Index, also known as `Pypi
<https://pypi.python.org/pypi>`__, in at the `Pypi PYDSM page
<https://pypi.python.org/pypi/pydsm>`__.

The pre-built documentation is also available in a dedicated `PyDSM
documentation <http://pythonhosted.org/pydsm/>`_ page on Pypi.

----

PyDSM is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your
option) any later version.

See file `COPYING` for further details.

Part of this code, limited to the `delsig` module, is ported from the
`DELSIG toolbox
<http://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox>`__
copyright by R. Schreier and licensed under the BSD license, as
specified in the corresponding files.

Distribution temporarily includes a patched version of the `CVXPY`
package by Tomas Tinoco de Rubira, that is currently discontinued,
being replaced by the `CVXPY` package by Steven Diamond and Eric Chu
and Stephen Boyd that provices more functionality and a different API.
This code is copyright by Tomas Tinoco de Rubira and licensed under the GPLv3+,
as specified in the corresponding files.

PyDSM is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY, without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.
