# MAGCPD

MAGCPD: A MATLAB code to calculate the Curie Point Depth and their 1D thermal modeling.

A MATLAB based codes for the estimation of the Curie Point Depth (CPD), associated with the depth to the bottom of magnetic source, by the inversion of magnetic anomalies using the method of Li et al. (2013) and defractal (Salem et al., 2014).

This code uses an algorithm based on the calculation of the top and the centroid of the magnetic source using the modified centroid approach due to 2D fractal magnetization model.

After the calculation of the CPD, an optional temperature at depth assuming a 1D steady state with heat generation could be calculate, introducing the thermal parameters required by the console. 
The functions requires as an input data:

a) GEOtiff image with the magnetic data or radially averaged power (or amplitude) spectrum (wavenumber (1/km), Log[P(k)] (or Log[A(k)])).

b) Flight altitude of the data acquisition in km.

Please watch the video tutorial for more details:
https://www.youtube.com/watch?v=Hgs5Yri0EDw

The synthetic data are:

1.- Magnetic data (synt_200.tif)

2.- Radially Averaged Amplitude Spectrum (RAAS.txt)

3.- Radially Averaged Power Spectrum (RAPS.txt)


This code is part of the paper: 
Carrillo-de la Cruz, J. L., Velázquez-Sánchez, S., Gómez-Rodríguez, D., Prol-Ledesma, R. M. MAGCPD: A MATLAB-based code to calculate the Curie point-depth involving the spectral analysis of aeromagnetic data. 

References:

Li, C.-F, Wang, J., Lin, J., Wang, T., 2013. Thermal evolution of the North Atlantic lithosphere: New constraints from magnetic anomaly inversion with a fractal magnetization model. Geochem. Geophys. Geosyst. 14 (12), 5078-5105. https://doi.org/10.1002/2013GC004896. 

Salem, A., Green, C., Ravat, D., Singh, K. H., East, P., Fairhead, J. D., Mogren, S., Biegert, E., 2014. Depth to Curie temperature across the central Red Sea from magnetic data using the defractal method. Tectonophysics 624–625, 75–86. https://doi.org/10.1016/j.tecto.2014.04.027. 
