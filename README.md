# MAGCPD

MAGCPD: MATLAB GUI to calculate the Curie Point Depth and their 1D thermal modeling.

A MATLAB based GUI for the estimation of the Curie Point Depth (CPD), associated with the depth to the bottom of magnetic source, by the inversion of magnetic anomalies using the modified centroid (Li et al. 2013) and defractal (Salem et al. 2014).

This code uses an algorithm based on the calculation of the top and the centroid of the magnetic source due to 2D fractal magnetization model.

After the calculation of the CPD, an optional temperature at depth profile is calculated assuming a 1D steady state with heat generation. For the temperature profile thermal parameters are required.

The GUI requires as an input data:

a) GEOtiff image with the magnetic data.
b) Flight altitude of the data acquisition in km.
c) Fractal parameters.
d) Selection interactive of the wavenumber range.
e) Thermal parameters for the temperature profiles

Data are attached and corresponds by:

1.- Synthetic Magnetic data (synt_200.tif)
2.- Example window (V9.tif)

The GUI can be installed as a standalone application using the . exe files. The process download the MATLAB runtime and install the application to used in your computer.

Also, the .mlapp are added to users that have license of MATLAB. This software is tested in MATLAB 2019a.


This code is part of the paper: 
Carrillo-de la Cruz, J. L., Velázquez-Sánchez, S., Gómez-Rodríguez, D., Prol-Ledesma, R. M. MAGCPD: A MATLAB-based GUI to calculate the Curie point-depth involving the spectral analysis of aeromagnetic data. Submit at Earth Science Informatics. 

References:

Li, C.-F, Wang, J., Lin, J., Wang, T., 2013. Thermal evolution of the North Atlantic lithosphere: New constraints from magnetic anomaly inversion with a fractal magnetization model. Geochem. Geophys. Geosyst. 14 (12), 5078-5105. https://doi.org/10.1002/2013GC004896. 

Salem, A., Green, C., Ravat, D., Singh, K. H., East, P., Fairhead, J. D., Mogren, S., Biegert, E., 2014. Depth to Curie temperature across the central Red Sea from magnetic data using the defractal method. Tectonophysics 624–625, 75–86. https://doi.org/10.1016/j.tecto.2014.04.027. 
