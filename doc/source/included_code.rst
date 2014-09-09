.. _included_code:

Code included in PyDSM
----------------------

Versions 0.9.x of PyDSM include the ``cvxpy`` package by Tomas Tinoco de
Rubira.

This is a discontinued library for modeling convex optimization
problems in Python. It has been superseeded by another library by
Steven Diamond, Eric Chu and Stephen Boyd that carries the same name,
but has a different API.

Since PyDSM was developed on the API of the former ``cvxpy``, in order
to support functionality the former ``cvxpy`` package is distributed
inside PyDSM as a temporary solution. In the long term, PyDSM will move
to use the newer ``cvxpy`` by Steven Diamond.

The ``cvxpy`` library provided with PyDSM was originally distributed
under the GPL3 license by its original author. See the source files
for further information. The PyDSM distribution installs it outside of
the ``pydsm`` package as ``cvxpy_tinoco`` in order to avoid name
conflicts with the newer ``cvxpy``.

The version shipped is the latest available from the original author
plus patches by Sergio Callegari and others. Some of these patches fix
minor errors with the code. A substantial set of patches is used to make
all the internal imports of the library relative, in order to make the
library relocatable to a new name.
