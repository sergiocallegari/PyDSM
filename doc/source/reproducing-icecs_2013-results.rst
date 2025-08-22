ICECS 2013 paper by Callegari
-----------------------------

This section illustrates how to replicate the results presented in the
paper [Cal13c]_.

To this aim, some sample code is provided in the directory
``Examples/ICECS-2013``.

The sample code reproduces the plots in Figs. 4 and 5 in the paper.
Furthermore, it outputs tabular data including that in Tables I.

With respect to the numeric data, note that the discrete time
simulations of the digital delta sigma modulator include some random
dithering. Thus the numeric data will never be identical to that in
the paper.

To start the examples, run:

``demo_dep_osr_hinf.py``
   to produce the plots in Fig. 4

``plot_dual_delsig.py``
   to produce the plots in Fig. 5 and data in the first part of table I

``plot_single_delsig.py``
   to produce plots similar to those in Fig. 5, but based on a reference
   case, and to produce the data in the second part of table I
