# PyDSM

A Python/Scipy toolbox for the design and simulation of ΔΣ Modulators, with emphasis on digital modulators.

Currently, the package focuses on tools for the design of the modulator NTF. It also includes a fast simulator for digital modulators. Furthermore, it includes a Python/Scipy port of some functions from the [DELSIG toolbox](http://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox) by R. Schreier.

The code includes an original NTF design technique fully described in the papers:

- Sergio Callegari, Federico Bizzarri *“Output Filter Aware Optimization of the Noise Shaping Properties of ΔΣ Modulators via Semi-Definite Programming”*, IEEE Transactions on Circuits and systems - Part I: Regular Papers, Vol. 60, N. 9, pp. 2352-2365. Sept. 2013. DOI: [10.1109/TCSI.2013.2239091](http://dx.doi.org/10.1109/TCSI.2013.2239091). Pre-print available on [ArXiv](http://arxiv.org/abs/1302.3020).

- Sergio Callegari, Federico Bizzarri *“Noise Weighting in the Design of ΔΣ Modulators (with a Psychoacoustic Coder as an Example),”* IEEE Transactions on Circuits and Systems - Part II: Express Briefs, Vol. 60, N. 11, pp. 756-760. Nov. 2013. DOI: [10.1109/TCSII.2013.2281892](http://dx.doi.org/10.1109/TCSII.2013.2281892). Pre-print available on [ArXiv](http://arxiv.org/abs/1309.6151).

If you find the code useful, *please cite these papers in your work*. In case you use the arxiv versions of the papers, make sure you cite the journal version and not the arxiv one.

## Code availability

The code is available on the Python Package Index, also known as [Pypi](https://pypi.python.org/pypi), at the [Pypi PYDSM page](https://pypi.python.org/pypi/pydsm).

No wheels are currently provided and the code builds from a source
distribution. On Linux, a blas library must currently be available at
the system level.

The code is routinely tested on Linux (AMD64) and occasionally on Windows 11 (again AMD64). It is also expected to work on MacOS, but currently untested.

## Documentation

The pre-built documentation is also available in a dedicated [PyDSM
documentation](http://pythonhosted.org/pydsm/) page on Pypi.

## Changelog

See [`CHANGELOG.rst`](doc/source/changelog.rst).

## Warning

Version 0.15.0 of PyDSM is substantially the same as version 0.14.0.0, with:

- modifications to the package build process, that now relies on modern Python practices (PEP 517);

- support for modern Python (>=3.10). Support for older Python is removed (even if there are still leftovers in the code itself);

- extension modules now built in Python 3 mode;

- a full list of dependencies provided in `pyproject.toml`, helping development with tools like `uv` or `pdm`.

Note that:

- Some warnings are to be expected running tests, with particular reference to the older `cvxpy` package by Tomas Tinoco de Rubiera that is currently vendored into PyDSM as `cvxpy_tdr`.

- Some tests may even fail, due to incorrect setup of the accuracy related parameters of the optimizers (may happen particularly with PICOS and SCS).

- The documentation now needs to be built using the `Makefile` in the `doc` directory. Only the html version of the documentation is expected to build. The documentation is substantially unchanged with respect to version 0.14.0.0, and so might be slightly misaligned with the current code. There may be also formatting issues due to the fact that the documentation source was originally written with a much older version of Sphinx in mind.

## Licensing information

PyDSM is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

See file [`COPYING`](COPYING) for further details.

Part of this code, limited to the `delsig` module, is ported from the
[DELSIG toolbox](http://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox) copyright by R. Schreier and licensed under the BSD license, as specified in the corresponding files.

Distribution temporarily includes a patched version of the `CVXPY` package by Tomas Tinoco de Rubira, that is currently discontinued, being replaced by the `CVXPY` package by Steven Diamond and Eric Chu and Stephen Boyd that provices more functionality and a different API. This code is copyright by Tomas Tinoco de Rubira and licensed under the GPLv3+, as specified in the corresponding files.

PyDSM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY, without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
