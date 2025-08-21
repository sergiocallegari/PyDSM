.. _included_code:

Code vendored in PyDSM
----------------------

Since version 0.9.0, PyDSM has vendored the original ``cvxpy`` package
by Tomas Tinoco de Rubira.

This was likely the first Python package to provide a modeling
framework for convex optimization problems. It has since been
discontinued and superseded by another package —also called ``cvxpy``—
developed by Steven Diamond, Eric Chu, and Stephen Boyd, which follows
a different API.

PyDSM was originally developed against Tinoco de Rubira’s ``cvxpy``.
When development of that package ceased and the new ``cvxpy`` emerged,
there arose a need to maintain compatibility and ensure a smooth
transition. The goals were:

* to preserve full reproducibility of published scientific results,
  avoiding significant variations caused by changes in the underlying
  optimization code;
* to allow meaningful comparison between the two codebases on the
  optimization problems relevant to PyDSM, both in terms of accuracy
  and computational efficiency.

For these reasons, PyDSM began supporting both versions of ``cvxpy``
by vendoring the original one. Tinoco de Rubira’s ``cvxpy`` was
originally distributed under the GPLv3 license. Its original codebase
has now disappeared from the web, but in PyDSM it is vendored as
:mod:`pydsm.cvxpy_tdr` to distinguish it from the newer ``cvxpy``. The
licensing and copyright information is fully preserved. The vendored
copy is based on the last release by the original author, with patches
by Sergio Callegari and others. While some patches fix minor issues,
most of them convert internal imports to relative imports, to enable
renaming and vendoring of the package without conflicts.

Now that the new ``cvxpy`` has matured, offers much better performance,
and alternatives such as ``PICOS`` have become available, the vendoring
of the original ``cvxpy`` is expected to be phased out in the near
future.
