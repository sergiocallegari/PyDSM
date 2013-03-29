The NTFdesign module
~~~~~~~~~~~~~~~~~~~~

.. module:: NTFdesign

Functions in this module are accessed by importing ``pydsm`` and
looking in the ``pydsm.NTFdesign`` module.

This module provides some functions for the design of the Noise
Transfer Function of ΔΣ modulators. There are both functions that are
specific to PyDSM and entry points to functions in :doc:`the delsig
module of PyDSM <delsig-module>`.

.. automodule:: pydsm.NTFdesign.delsig
.. module:: NTFdesign.delsig

.. autofunction:: pydsm.NTFdesign.delsig.synthesizeNTF
.. autofunction:: pydsm.NTFdesign.delsig.synthesizeChebyshevNTF

.. automodule:: pydsm.NTFdesign.filter_based
.. module:: NTFdesign.filter_based

**Key functions**

.. autofunction:: pydsm.NTFdesign.filter_based.synthesize_ntf_from_filter_imp

**Auxiliary Functions**

.. autofunction:: pydsm.NTFdesign.filter_based.synthesize_ntf_from_q0
.. autofunction:: pydsm.NTFdesign.filter_based.q0_from_filter_imp_response
.. autofunction:: pydsm.NTFdesign.filter_based.quantization_noise_gain
