# GPR data processing
Processing of ground penetrating radar (GPR) data using Matlab and the CREWES library, an open-source (mostly) seismic data processing library. Amongst other things, this is an implementation of the method described in [Harlan et al., Geophysics, vol 49, no. 11, 1984](https://library.seg.org/doi/pdf/10.1190/1.1441600). It is applied to georadar data to filter horizontal beds and diffractions and to estimate wave velocities.

# Work in progress

Decembre 10 2017: The repo is a total mess. It is not quite usefull for anyone as is. This was a quick and dirty implementation of Harlan et al.'s method I did in 2015. It seemed to work pretty well on good GPR data (see next section) but is totally unoptimized.

1) Cleaning up the code documentation and commenting. A lot of m-files are not relevent to the implementation either.
2) I had to modify the code of the fk-migration in the CREWES library for the migration to work on time series not starting at t=0 but I did not save this modification at the time -_-. Need to figure it out anew.
3) The code makes use of the kernel density estimation (KDE) for some deconvolution. I used code by [Z. I. Botev](https://arxiv.org/abs/1011.2602). As I recall, this KDE step was the bottle neck of the implementation in terms of running time. This has to be investigated for optimization

# Background 

This was a quick and dirty implementation of Harlan et al.'s method that seemed to work on GPR data. Here is the GPR profile I used for tests.

<img src="graphique_focalisation.tif" width="900px"/>

![alt text](https://github.com/dabana/gpr_processing/blob/master/graphique_focalisation.tif)

