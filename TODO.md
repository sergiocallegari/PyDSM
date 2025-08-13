# To Do list for the pydsm code-base

## Perspective

### Framework

- [ ] Set up continuous CI/CD
- [ ] Set up things to have signed release TAGS

### Options for modelers and optimizers

- [ ] Improve option handling for modelers and optimizers
  - Consider moving to an object oriented interface?

### Performance

- [ ] Consider changing default modeler/optimizer
  - Maybe cvxpy/scs
- [ ] Consider Clarabel and its performance

### DELSIG

- [ ] Maybe have a plain port of delsig rather than implementation of part of its functions on top of pydsm optimizers and then propose some delsig "emulation" within pydsm
- [ ] Look into quadrature modulators
- [ ] Implement `findPattern`
- [ ] Implement `calculateSNR`
- [ ] Implement `simulateSNR`
- [ ] Consider the following suggestion for the CLANS method, that was submitted as a patch to the clans Matlab code: for better convergence, transform the errors in the minimization functions (`dsclansObj6a`, `dsclansObj6b`) into quadratic errors. Steps for the original MATLAB code
  1. Open `clans6.m`;
  2. in sub-function `dsclansObj6a`, after `f = abs(evalTF(H,exp(1i*pi/OSR)));` type `f=f*f`;
  3. in sub-function `dsclansObj6a`, after `g = sum(abs(impulse(H,100))) -1 - Q;` type `g=g*g`.
      In this way the target value for `Q` will likely be reached quite accurately at low OSR.

### SDTOOLBOX

- [ ] Implement `calcSNR`
- [ ] Implement `sinusx`

### Simulator

- [ ] If code based on external cblas is retained, reuse common code between `_simulateDSM_scipy_blas` and `_simulateDSM_cblas`
- [ ] Consider other blas or blas like options for the simulator (accelerated; fblas; different matrix representations (C, fortran); or different approaches for the simulator)

### DSM as an euristic optimizer

- [ ] Implement quadratic form
- [ ] Implement scanner with quadratic form


### Documentation

- [ ] Provide a tutorial
- [ ] Provide notebook version of the examples
- [ ] Maybe use the versionadded/versionchanged/deprecated sphinx features to document API changes
- [ ] In the examples, the numeric results in the ICECS 2012 paper are not consistent with what the sample script returns. Specifically, the SNR values evaluated via the NTF on the paper are better than those obtained by the code by a factor 2. In the original version of the code, the results where identical due to a missing factor 2 in the quantization noise gain evaluation. Probably, this error slipped in because it makes the NTF-based results and the time-domain-simulation results more similar. In fact, the *real* modulator behavior is better than predicted by the NTF model because in the specific example the quantization noise is not white, but for some reason slightly blue (more power at higher frequencies). A note should probably be added to the code and to the ArXiv version of the ICECS paper.

### Refactoring

- [ ] Consider moving from plain `zpk` and `ba` representations to named tuple model or to LTI for filters


## Deferred

### Framework

- [x] Avoid local links in `README.md` that may fail from `pypi`.
  - A `README-PyPI.md` file is now generated dynamically from `README.md` in `setup.py`. This is not very canonical, but appears to work.
- [ ] Make some modelers optional
- [ ] Provide interface to expose available features

### Legacy

- [ ] Avoid `__import__` statements and follow https://docs.scipy.org/doc/scipy/reference/
- [ ] Completely remove vendored `cvxpy` package by Tomas Tinoco de Rubiera
  - Superseded by `cvxpy` and `PICOS`
  - Relies on `np.matrix` that throws warnings and is deprecated
  - Needs patches anyway with recent NumPy due to use of `np.NaN`
- [ ] Consider removal of `_simulateDSM_cblas`
  - Means relying only on the cblas provided by numpy/scipy
  - Simplifies building, but currently the scipy cblas is slightly less efficient on Linux
  - A notable issue to consider is the lack of `cblas.h` in the openblas provided by Numpy/Scipy
- [ ] Remove cruft meant to support python 2.7
  - Specifically, get rid of `six`
- [ ] Consider pyupgrade for modernizing code
- [ ] Consider ruff for checking the code
- [ ] Get rid of deprecated `matplotlib.mlab` in the ICECS-2013 examples
- [ ] change `logspace` into `geomspace`
- [ ] avoid inline function definitions

### Accuracy

- [ ] A recent code change has lead to removal of a matrix inversion noticing that the matrix is orthonormal, and that transposition should suffice.  However this seems to cause come accuracy issues on some platforms (AMD?). Check if this still the case.

### Simulator

- [ ] Consider forcing some type casts in the cython based simulator code (e.g. `np.double`)
- [ ] Consider using `numpy.require` to assure that arrays have the right properties for interfacing with cython code


## For the 0.15.0 release

Provide minimal changes to support deployment on modern Python

### Framework

- [x] Move to a modern build system (no direct invocation of `setup.py`)
  - Keep using `setuptools`
  - Build using `uv build` or maybe `pdm build`
- [x] Use proper lock files
  - OK with `uv`, `pdm`
- [x] Move code to `src`
- [x] Avoid moving tests to `tests`
  - since tests should be part of the pakage (invokable via a `pydsm.test()` function)
- [x] Update support for benchmarking
  - benchmarking should also be a part of the package (invokable via a `pydsm.bench()` function)
- [x] Switch from `nose` to `unittest` or `pytest`
  - Adpoted `pytest` for alignment with `numpy`.
- [x] Ensure that `cython` is a *build* requirement only
- [x] Ensure that the `sphynx` documentation builds
  - Some rendering issues remain, only the HTML documentation is currently tested
- [x] Support latest numpy/scipy
  - Warnings remain due to `tinoco_cvxpy`
- [x] Support latest python
  - Currently 3.13.x
- [x] get rid of `importlib_resources`
  - The functionality is directly provided in recent Python
- [x] Clean up uv/setuptools build warnings
    - [x] Use SPDX license expression for `project.license` or `project.license-files` in place of classifiers. See [license Instructions for `pyproject.toml`](https://packaging.python.org/en/latest/guides/writing-pyproject-toml/#license) for details. Adopt `License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)`
	- [x] Fix inclusion of benchmark material `pydsm.NTFdesign.benchmarks` and `pydsm.delsig.benchmarks.Data` adding it to the `packages` configuration field
	- [x] Fix inclusion of test material `pydsm.delsig.tests.Data`
- [x] Avoid inclusion of unneeded files in the sdist (template, docs, examples)
    - Obtained automatically by moving the source package into `src`
- [x] Use `setuptools_scm` rather than custom solution


### Code examples

- [x] Port the code examples to proper Python 3
