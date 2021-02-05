# MAGCPD

MAGCPD: MATLAB GUI to calculate the Curie Point Depth and their 1D thermal modeling.

A MATLAB based GUI for the estimation of the Curie Point Depth (CPD), associated with the depth to the bottom of magnetic source, by the inversion of magnetic anomalies using the modified centroid (Li et al. 2013) or the defractal (Salem et al. 2014) methods.

This code uses an algorithm based on the calculation of the top and the centroid of the magnetic source due to 2D fractal magnetization model.

After the calculation of the CPD, an optional temperature at depth profile is calculated assuming a 1D steady state with heat generation. For the temperature profile thermal parameters are required.

The GUI requires as an input data:

a) GEOtiff image with the magnetic data. Also, old Geosoft 2-byte (signed integer) binary grid are other format to the program. It is important to notice that the square window need to be exported from Oasis Montaj as “Geosoft short” grid type.

b) Flight altitude of the data acquisition in km.

c) Fractal parameters.

d) Interactive selection of the wavenumber range.

e) Thermal parameters for the temperature profiles.


Data are attached and corresponds by:

1.- Synthetic Magnetic data (synt_200.tif).

2.- Example window (V9.tif).

3.- Example for .grd files (example.grd). There is two images (df_example.png and mc_example.png) showing the wavenumber range and fractal parameter used for this example

The GUI can be installed as a standalone application using the . exe files. The process download the MATLAB runtime and install the application to used in your computer.

Also, the .mlapp are added to users that have license of MATLAB. This software is tested in MATLAB 2020a. ATTENTION: selectdata.m need to be in the same folder of the .mlapp applications for matlab running.


This code is part of the paper: 
Carrillo-de la Cruz J. L.,  Prol-Ledesma R. M., Velázquez-Sánchez S., Gómez-Rodríguez D. (2020).  MAGCPD: A MATLAB-based GUI to calculate the Curie point-depth involving the spectral analysis of aeromagnetic data. Earth Science Informatics. https://doi.org/10.1007/s12145-020-00525-x

References:

Li, C.-F, Wang, J., Lin, J., Wang, T., 2013. Thermal evolution of the North Atlantic lithosphere: New constraints from magnetic anomaly inversion with a fractal magnetization model. Geochem. Geophys. Geosyst. 14 (12), 5078-5105. https://doi.org/10.1002/2013GC004896. 

Salem, A., Green, C., Ravat, D., Singh, K. H., East, P., Fairhead, J. D., Mogren, S., Biegert, E., 2014. Depth to Curie temperature across the central Red Sea from magnetic data using the defractal method. Tectonophysics 624–625, 75–86. https://doi.org/10.1016/j.tecto.2014.04.027.

## Important note: after the instalation of the .exe file, the first time to run the program it may take some minutes, please be patient! 

## Important change for the new commit (02/10/2020):
Recently, a researcher contact us to show us an error in the equation 8 of the paper, this due to an use of an old version of the paper "Martos YM, Catalán M, Galindo‐Zaldivar J (2019) Curie Depth, heat ﬂux, and thermal subsidence reveal the Paciﬁc mantle outﬂow through the Scotia Sea. J Geophys Res: Solid Earth 124: 10,735–10,751. https://doi.org/10.1029/2019JB017677". The equation is: 

<img src="https://render.githubusercontent.com/render/math?math=\Delta Z_{b} =\sqrt{ 2 \cdot \Delta Z_{0}^2 %2B \Delta Z_{t}^2 }">.

And the correct equation:

<img src="https://render.githubusercontent.com/render/math?math=\Delta Z_{b} = \sqrt{ (2 \Delta Z_{0}) ^2 %2B \Delta Z_{t}^2 }">.

Our program is updated with the correct equation

# Version 2.0 is ready

This new version has new improvements in the calculation of the radially averaged amplitude spectrum.  Also, the magnetic data grid resolution needs to be specified in km.
