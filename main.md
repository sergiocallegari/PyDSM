# Python/Scipy tools for the design and simulation of ΔΣ modulators #

This is a Python/Scipy toolbox for the design and simulation of Digital Delta Sigma Modulators.

Currently, it focuses on tools for the design of the modulator NTF.  It also includes a fast simulator for digital modulators. Furthermore, it includes a Python/Scipy port of some functions from the DELSIG toolbox by R. Schreier.

As a highlight, the code includes an original NTF design technique fully described in the papers:

  * Sergio Callegari, Federico Bizzarri _“Output Filter Aware Optimization of the Noise Shaping Properties of ΔΣ Modulators via Semi-Definite Programming”_, IEEE Transactions on Circuits and Systems - Part I: Regular Papers, Vol. 60, N. 9, pp. 2352-2365. Sept. 2013. DOI: [10.1109/TCSI.2013.2239091](http://dx.doi.org/10.1109/TCSI.2013.2239091). Pre-print available on [ArXiv](http://arxiv.org/abs/1302.3020).

  * Sergio Callegari, Federico Bizzarri _“Noise Weighting in the Design of ΔΣ Modulators (with a Psychoacoustic Coder as an Example),”_ IEEE Transactions on Circuits and Systems - Part II: Express Briefs, Vol. 60, N. 11, pp. 756-760. Nov. 2013. DOI: [10.1109/TCSII.2013.2281892](http://dx.doi.org/10.1109/TCSII.2013.2281892). Pre-print available on [ArXiv](http://arxiv.org/abs/1309.6151).

**If you find the code useful, please, cite the papers in your work**

The code base also includes sample code to replicate the results in the above papers and in:

  * Sergio Callegari, Federico Bizzarri _“Should ΔΣ modulators used in AC motor  drives be adapted to the mechanical load of the motor?”_, Proceedings of the 19th IEEE International Conference on Electronics, Circuits and Systems (ICECS), 2012, pp. 849 - 852. DOI: [10.1109/ICECS.2012.6463619](http://dx.doi.org/10.1109/ICECS.2012.6463619). Pre-print available on [arXiv](http://arxiv.org/abs/1302.7172).

  * Sergio Callegari, _“Should ΔΣ modulators used in AC motor drives be adapted to the mechanical load of the motor?”_, Proceedings of the 20th IEEE International Conference on Electronics, Circuits and Systems (ICECS), 2013, pp. 589 - 592. DOI: [10.1109/ICECS.2013.6815483](http://dx.doi.org/10.1109/ICECS.2013.6815483).

While downloading the software, note that the wiki on this site contains an [introduction](Introduction.md) to the project and some automatically generated [toolbox documentation](http://wiki.pydsm.googlecode.com/git/Documentation/index.html) for further information. Why not taking a look at it? With this, you'll be immediately ready to use the code.

## New! Version 0.10.1 is here! ##

This version includes the following changes:

  * Fixes to build on MacOs
  * Documentation instructions go get started on MacOs
  * Some documentation fixes

## Licensing notice ##

Code License is GPL V3. Up to version 0.8.1 code licensing was indicated as BSD, but this forbade distribution with the optimization toolboxes that are used from PyDSM.