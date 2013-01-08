How to replicate the results in the TCAS-I paper
------------------------------------------------

This section illustrates how to replicate the results presented in the
paper

    Sergio Callegari, Federico Bizzarri *"Output Filter Aware
    Optimization of the Noise Shaping Properties of ΔΣ Modulators via
    Semi-Definite Programming"*, IEEE Transactions on Circuits and
    Systems - Part I: Regular Papers. To appear in 2013.

To this aim, some sample code is provided in the directory
``Examples/TCAS1-2013``.

The research paper provides 3 examples, in Sections V-a, V-b and V-c.

Example in section V-a
''''''''''''''''''''''

This example refers to the design of a Digital ΔΣ modulator that is
followed by a low pass filter in charge of removing the quantization
noise.

The corresponding code is provided in the file ``demo-lp.py``. Once
PyDSM and all its pre-requisites are installed, this can be started
directly by opening a shell (command prompt) and typing::

  python demo-lp.py

Alternatively, the script can be opened in Spyder and launched from
there.

The code runs, showing some intermediate output from the
optimizer. Then it provides the graphical output that is reported in
Figures 7 and 9 in the TCAS-I paper.

The code also provides a further graph, obtained by the time domain
simulation of the modulator. This includes 3 curves:

#. The output of a ΔΣ modulator passed through the filter, where the
   modulator is designed by the technique in the TCAS-I paper.

#. The output of a ΔΣ modulator passed through the filter, where the
   modulator is designed by the ``synthesizeNTF`` function from the
   DELSIG toolbox.

#. The input of the modulator passed through the filter.

This latter plot is meant to provide a graphical illustration that the
modulator obtained by the proposed design methodology is actually
working correctly. All the three plots should overlap almost
perfectly. However, zooming in, little differences should become
apparent, with the curve from the optimized modulator following the
curve given by the input signal slightly better.

The code also provides some numerical output, corresponding to the SNR
values reported in the research paper, obtained both with the
linearized model and with the actual time-domain simulation of the
modulator.

Note that the code does not provide any equivalent of Figure 8. For
this, it is necessary to introduce a loop so that the NTF is optimized
for different modulator orders. This is easy to achieve.

Neither the current code provides an equivalent of Figure 10. This
is anyway easy to achieve by using the ``plotPZ`` function in the
``pydsm.delsig`` module.

As a final remark, the proposed example code is not particularly
elegant, but should be rather easy to read, also thanks to the many
comments.


Example in section V-b
''''''''''''''''''''''

This example refers to the design of a Digital ΔΣ modulator that is
followed by a bandpass filter in charge of removing the quantization
noise.

The corresponding code is provided in the file ``demo-bp.py``. Once
PyDSM and all its pre-requisites are installed, this can be started
directly by opening a shell (command prompt) and typing::

  python demo-bp.py

Alternatively, the script can be opened in Spyder and launched from
there.

The code runs, showing some intermediate output from the
optimizer. Then it provides the graphical output that is reported in
Figures 11 and 13 in the TCAS-I paper.

The code also provides a further graph, obtained by the time domain
simulation of the modulator. This includes 3 curves:

#. The output of a ΔΣ modulator passed through the filter, where the
   modulator is designed by the technique in the TCAS-I paper.

#. The output of a ΔΣ modulator passed through the filter, where the
   modulator is designed by the ``synthesizeNTF`` function from the
   DELSIG toolbox.

#. The input of the modulator passed through the filter.

This latter plot is meant to provide a graphical illustration that the
modulator obtained by the proposed design methodology is actually
working correctly. All the three plots should overlap almost
perfectly. However, zooming in, little differences should become
apparent, with the curve from the optimized modulator following the
curve given by the input signal slightly better.

The code also provides some numerical output, corresponding to the SNR
values reported in the research paper, obtained both with the
linearized model and with the actual time-domain simulation of the
modulator.

Note that the code does not provide any equivalent of Figure 12. For
this, it is necessary to introduce a loop so that the NTF is optimized
for different modulator orders. This is easy to achieve.

Nor the current code provides an equivalent of Figure 10. This
is anyway easy to achieve by using the ``plotPZ`` function in the
``pydsm.delsig`` module.

As a final remark, the proposed example code is not particularly
elegant, but should be rather easy to read, also thanks to the many
comments.

**Warning** The optimization of the modulator NTF can be particularly
time consuming in Windows, where only a 32 bit computation environment
is available and where some math libraries are not fully optimized to
the specific hardware platform in CVXOPT. For maximum performance, use
a modern 64 bit Linux distribution on a multi core CPU with fast
floating point math and instructions set with SIMD (SSEx) extensions
and possibly Advanced Vector Extensions (AVX).


Example in section V-c
''''''''''''''''''''''

This example refers to the design of a Digital ΔΣ modulator that is
followed by a two band filter in charge of removing the
quantization noise.

The corresponding code is provided in the file
``demo-multiband.py``. Once PyDSM and all its pre-requisites are
installed, this can be started directly by opening a shell (command
prompt) and typing::

  python demo-multiband.py

Alternatively, the script can be opened in Spyder and launched from
there.

The code runs, showing some intermediate output from the
optimizer. Then it provides the graphical output that is reported in
Figures 15 and 17 in the TCAS-I paper.

Note that the sample code provides no comparison to Nagahara's min-max
design strategy.

The code also provides a further graph, obtained by the time domain
simulation of the modulator. This includes 2 curves:

#. The output of a ΔΣ modulator passed through the filter, where the
   modulator is designed by the technique in the TCAS-I paper.

#. The input of the modulator passed through the filter.

This latter plot is meant to provide a graphical illustration that the
modulator obtained by the proposed design methodology is actually
working correctly. The two plots should overlap almost
perfectly.

The code also provides some numerical output, corresponding to the SNR
values reported in the research paper, obtained both with the
linearized model and with the actual time-domain simulation of the
modulator.

Note that the code does not provide any equivalent of Figure 16. For
this, it is necessary to introduce a loop so that the NTF is optimized
for different modulator orders. This is easy to achieve.

Nor the current code provides an equivalent of Figure 18. This
is anyway easy to achieve by using the ``plotPZ`` function in the
``pydsm.delsig`` module.

As a final remark, the proposed example code is not particularly
elegant, but should be rather easy to read, also thanks to the many
comments.

**Warning** The optimization of the modulator NTF can be particularly
time consuming in Windows, where only a 32 bit computation environment
is available and where some math libraries are not fully optimized to
the specific hardware platform in CVXOPT. For maximum performance, use
a modern 64 bit Linux distribution on a CPU with fast floating point
math and instructions set with SIMD (SSEx) extensions and possibly
Advanced Vector Extensions (AVX).
