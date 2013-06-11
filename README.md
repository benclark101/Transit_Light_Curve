# Transit Light Curve Analysis Pipeline
Allows the creation and analysis of transiting exoplanet light curves directly from FITS files

# This routine can:
* Perform aperture photometry on FITS files
* Create basic light curves using calibration stars (and a super-calibrator)
* Calculate the radius, period, inclination angle, semi-major axis as well as the density and mass of the  exoplanet (providing stellar information is known and MpSin(i) from radial velocity measurements)
* Fit advanced models to the light curve and obtain transit midpoint time and luminosity dip
* Perform Monte-Carlo Mandel-Agol limb darkening model fitting 
* Create an Observed-Calculated plot for the exoplanet (currently has to be done manually)

# Requirements:
Python 2.7
The following libraries:
* [pyFITS](http://www.stsci.edu/institute/software_hardware/pyfits)
* [NumPy](http://www.numpy.org/)
* csv
* os
* math
* [matplotlib](http://matplotlib.org/)
* [SciPy](http://www.scipy.org/)
* random
I recommend you download [PythonXY](https://code.google.com/p/pythonxy/)

# To install:
Download repository as a zip and extract it to a folder of your choice

# To Run:
1. Edit "userinput.py" and enter the location of the FITS file you wish to analyse
2. Enter host star properties in the same file
3. Run "main.py" and follow instructions when prompted

