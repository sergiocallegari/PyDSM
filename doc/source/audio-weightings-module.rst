The audio-weightings module
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Functions in this module are accessed by importing ``pydsm`` and looking in the
``pydsm.audio_weightings`` module.

.. automodule:: pydsm.audio_weightings

Key functions
.............

.. autofunction:: pydsm.audio_weightings.a_weighting
.. autofunction:: pydsm.audio_weightings.b_weighting
.. autofunction:: pydsm.audio_weightings.c_weighting
.. autofunction:: pydsm.audio_weightings.d_weighting
.. autofunction:: pydsm.audio_weightings.f_weighting

Data elements
.............

* ``a_zpk`` - filter in zpk form implementing the A-weighting
* ``b_zpk`` - filter in zpk form implementing the B-weighting
* ``c_zpk`` - filter in zpk form implementing the C-weighting
* ``d_zpk`` - filter in zpk form implementing the D-weighting
* ``f_zpk`` - filter in zpk form implementing the F-weighting

These filters are unnormalized, meaning that the gain at 1 kHz may be arbitrary.
