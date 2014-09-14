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
.. autofunction:: pydsm.NTFdesign.delsig.clans

.. automodule:: pydsm.NTFdesign.weighting
.. module:: NTFdesign.weighting

**Key functions**

.. autofunction:: pydsm.NTFdesign.weighting.synthesize_ntf_from_noise_weighting

**Auxiliary functions**

.. autofunction:: pydsm.NTFdesign.weighting.q0_from_noise_weighting
.. autofunction:: pydsm.NTFdesign.weighting.synthesize_ntf_from_q0


.. automodule:: pydsm.NTFdesign.filter_based
.. module:: NTFdesign.filter_based

**Key functions**

.. autofunction:: pydsm.NTFdesign.filter_based.synthesize_ntf_from_filter

**Auxiliary functions**

.. autofunction:: pydsm.NTFdesign.filter_based.q0_from_filter
.. autofunction:: pydsm.NTFdesign.filter_based.quantization_noise_gain

.. automodule:: pydsm.NTFdesign.psychoacoustic
.. module:: NTFdesign.psichoacoustic

**Key functions**

.. autofunction:: pydsm.NTFdesign.psychoacoustic.synthesize_ntf_dunn
.. autofunction:: pydsm.NTFdesign.psychoacoustic.synthesize_ntf_from_audio_weighting

**Auxiliary functions**

.. autofunction:: pydsm.NTFdesign.psychoacoustic.dunn_optzeros
.. autofunction:: pydsm.NTFdesign.psychoacoustic.dunn_optzeros_cplx

.. automodule:: pydsm.NTFdesign.minmax
.. module NTFdesign.minmax

.. autofunction:: pydsm.NTFdesign.minmax.synthesize_ntf_minmax
