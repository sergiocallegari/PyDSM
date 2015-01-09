Reference guide
---------------

The current codebase comprises both code that is specific to PyDSM and
code ported from the DELSIG_ Toolbox by R. Schreier. This is reflected
in the module organization of PyDSM.

Currently, the DELSIG port is rather limited and just includes what is
useful for comparison to the specific PyDSM algorithms.

Right now, PyDSM is mostly specialized in the design of Noise Transfer
Functions for Digital ΔΣ Modulators and in the simulation of the
resulting modulators.

.. toctree::
   :maxdepth: 1

   delsig-module
   ntf-design-module
   simulation-module
   audio-weightings-module
   iso226-module

.. toctree::
   :maxdepth: 2

   utility-modules

Error handling
..............

.. toctree::
   :maxdepth: 1

   pydsm.exceptions



.. _DELSIG:
   http://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox
