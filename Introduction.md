# Introduction #

PyDSM is a Python Delta Sigma Modulator toolbox. It contains tools for experimenting with ΔΣ modulators. Currently, it is still relatively reduced and focuses mostly on the experimentation of different techniques for the design of the modulator Noise Transfer Function (NTF). Furthermore, the current codebase contains means to simulate a generic digital modulator.

PyDSM is under development and shall be enriched with further functionalities in a near future.

Highlights of the toolbox consist in:

  * Some routines ported from the very well known [DELSIG Toolbox](http://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox) for Matlab by R. Schreier

  * The NTF design algorithm presented in:

  * Sergio Callegari, Federico Bizzarri _“Output Filter Aware Optimization of the Noise Shaping Properties of ΔΣ Modulators via Semi-Definite Programming”_, IEEE Transactions on Circuits and Systems - Part I: Regular Papers, Vol. 60, N. 9, pp. 2352-2365. Sept. 2013. DOI: [10.1109/TCSI.2013.2239091](http://dx.doi.org/10.1109/TCSI.2013.2239091). Pre-print available on [ArXiv](http://arxiv.org/abs/1302.3020).

  * Sergio Callegari, Federico Bizzarri _“Noise Weighting in the Design of ΔΣ Modulators (with a Psychoacoustic Coder as an Example),”_ IEEE Transactions on Circuits and Systems - Part II: Express Briefs, Vol. 60, N. 11, pp. 756-760. Nov. 2013. DOI: [10.1109/TCSII.2013.2281892](http://dx.doi.org/10.1109/TCSII.2013.2281892). Pre-print available on [ArXiv](http://arxiv.org/abs/1309.6151).

  * And more...

If you find this code useful, _please consider citing these papers in your work._

See the automatically generated [toolbox documentation](http://pythonhosted.org//pydsm) for further information.