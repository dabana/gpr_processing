# GPR data processing
Processing of ground penetrating radar (GPR) data using Matlab and the CREWES library, an open-source (mostly) seismic data processing library. Amongst other things, this is an implementation of the method described in [Harlan et al., Geophysics, vol 49, no. 11, 1984](https://library.seg.org/doi/pdf/10.1190/1.1441600). It is applied to georadar data to filter horizontal beds and diffractions and to estimate wave velocities.

# Work in progress

The repo is a total mess. This was a quick and dirty implementation of Harlan et al.'s method that seemed to work on GPR data (see next section) but is totally unoptimized and not quite 

1) Cleaning up the code documentation and commenting
2) I had to modify the code of the fk-migration in the CREWES library for the migration to work on time series not starting at t=0 but I did not save this modification at the time -_-. Need to figure it out anew.
3) The code makes use of the Kernel Deconvolution 

# Background 

The repository was created in december 2015. This was a quick and dirty implementation of Harlan et al.'s method that seemed to work on GPR data. Here is 

