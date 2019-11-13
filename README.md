# MAGCPD

MAGCPD: A MATLAB functions to calculate the Curie Point Depth and their 1D thermal modeling.

A MATLAB based codes for the estimation of the Curie Point Depth (CPD), associated with the depth to the bottom of magnetic source, by the inversion of magnetic anomalies.

This code uses an algorithm based on the calculation of the top and the centroid of the magnetic source using the modified centroid approach due to 2D fractal magnetization model.

After the calculation of the CPD, an optional temperature at depth assuming a 1D steady state with heat generation could be calculate, introducing the thermal parameters required by the console. 
The functions requires as an input data:

a) GEOtiff image with the magnetic data or radially averaged power (or amplitude) spectrum (wavenumber (1/km), Log[P(k)] (or Log[A(k)])).

b) Flight altitude of the data acquisition in km.

