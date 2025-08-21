TCAS-II 2013 paper by Callegari and Bizzarri
--------------------------------------------

This section illustrates how to replicate the results presented in the
paper [Cal13b]_.

To this aim, some sample code is provided in the directory
``Examples/TCAS2-2013``.

The research paper provides 4 examples, in Sections IV-a to V-e.

Example in section IV-a
'''''''''''''''''''''''

This example refers to the design of a Digital ΔΣ modulator for low
pass signals with a mere specification of the signal bandwith (which
is specified as a fraction of the modulator sample frequency via the
oversampling ratio).

The corresponding code is provided in the file ``demo_lp_brickwall.py``. Once
PyDSM and all its pre-requisites are installed, this can be started
directly by opening a shell (command prompt) and typing::

  python demo_lp_brickwall.py

Alternatively, the script can be opened in Spyder and launched from
there.

The code runs, showing some intermediate output from the
optimizer. Then it provides the graphical output that is delivered in
Figure 2 in the TCAS-II paper.

The proposed example code is not particularly elegant, but should be
rather easy to read, also thanks to the many comments.


Example in section IV-b
'''''''''''''''''''''''

This example refers to the design of a multiband Digital ΔΣ modulator.

The corresponding code is provided in the file
``demo_multiband_brickwall.py``. Once PyDSM and all its pre-requisites
are installed, this can be started directly by opening a shell
(command prompt) and typing::

  python demo_multiband_brickwall.py

Alternatively, the script can be opened in Spyder and launched from
there.

The code runs, showing some intermediate output from the
optimizer. Then it provides the graphical output that is delivered in
Figure 3 in the TCAS-II paper.

The proposed example code is not particularly
elegant, but should be rather easy to read, also thanks to the many
comments.

Example in section IV-c
'''''''''''''''''''''''

This example refers to the design of a Digital ΔΣ modulator that is
followed by a non ideal filter in charge of removing the quantization
noise.

The corresponding code is provided in the file
``demo_lp_filter.py``. Once PyDSM and all its pre-requisites are
installed, this can be started directly by opening a shell (command
prompt) and typing::

  python demo_lp_filter.py

Alternatively, the script can be opened in Spyder and launched from
there.

The code runs, showing some intermediate output from the
optimizer. Then it provides the graphical output that is delivered in
Figure 5 in the TCAS-II paper.

The proposed example code is not particularly elegant, but should be
rather easy to read, also thanks to the many comments.

Example in section IV-d
'''''''''''''''''''''''

This example refers to the design of a Digital ΔΣ modulator for low-pass
signal that is capable of delivering a reduced quantization noise close to
dc.

The corresponding code is provided in the file
``demo_lp_low_dc_noise.py``. Once PyDSM and all its pre-requisites are
installed, this can be started directly by opening a shell (command
prompt) and typing::

  python demo_lp_low_dc_noise.py

Alternatively, the script can be opened in Spyder and launched from
there.

The code runs, showing some intermediate output from the
optimizer. Then it provides the graphical output that is delivered in
Figure 6 in the TCAS-II paper.

The proposed example code is not particularly elegant, but should be
rather easy to read, also thanks to the many comments.

Example in section IV-e
'''''''''''''''''''''''

This example refers to the design of a Digital ΔΣ modulator for audio
signals capable of shaping the in-band quantization noise according to
a psychoacoustic weighting function, so that it is minimally audible.

The corresponding code is provided in the file
``demo_psychoacoustic.py``. Once PyDSM and all its pre-requisites are
installed, this can be started directly by opening a shell (command
prompt) and typing::

  python demo_psychoacoustic.py

Alternatively, the script can be opened in Spyder and launched from
there.

The code runs, showing some intermediate output from the
optimizer. Then it provides the graphical output that is delivered in
Figure 7 in the TCAS-II paper.

The proposed example code is not particularly elegant, but should be
rather easy to read, also thanks to the many comments.
