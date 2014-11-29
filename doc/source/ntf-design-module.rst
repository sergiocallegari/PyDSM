The NTFdesign module
~~~~~~~~~~~~~~~~~~~~

.. automodule:: pydsm.NTFdesign

   Functions in this module are accessed by importing ``pydsm`` and
   looking in the ``pydsm.NTFdesign`` module.

   This module provides some functions for the design of the Noise
   Transfer Function of ΔΣ modulators. There are both functions that are
   specific to PyDSM and entry points to functions in :doc:`the delsig
   module of PyDSM <delsig-module>`.

   The module includes the following submodules

   .. toctree::
      :maxdepth: 1

      ntf-design-weighting
      ntf-design-minmax
      ntf-design-delsig
      ntf-design-psychoacoustic
      ntf-design-legacy
      ntf-design-filter-based


   Key functions promoted at module top level
   ..........................................

   .. autofunction:: ntf_schreier
   .. autofunction:: ntf_chebyshev
   .. autofunction:: ntf_clans
   .. autofunction:: ntf_fir_weighting
   .. autofunction:: ntf_fir_minmax
   .. autofunction:: ntf_dunn
   .. autofunction:: ntf_fir_audio_weighting
   .. autofunction:: quantization_noise_gain
