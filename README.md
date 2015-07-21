# HG1G2 tools

This is a Python package for fitting the H,G1,G2 magnitude system to
asteroid photometric data.

This package is still under heavy developement and anything is
susceptible to change until further notice.

## Installation

Running `python setup.py install` should be enough. NumPy is required.

## Usage
### Example

```
import HG1G2tools

data = loadtxt("44_Nysa.dat", skiprows=1)
fit, covmatrix = HG1G2tools.fit_HG1G2(data)
error_estimates = HG1G2tools.simulate_errors(fit, covmatrix)

```


### Fitting H, G1, G2

The function `fit_HG1G2` performs the H, G1, G2 fit. The default form is

`fit, covmatrix = fit_HG1G2(data, weights=None, degrees=True)`

The parameter `data` needs to be an `(N,2)` shaped NumPy array, where `N` is the number
of data points. The first column should have the phase angle and the second column the
observed magnitude.

The optional parameter `degrees` is a boolean, (default `True`), specifying whether the
phase angle in `data` is given in degrees or radians.

The optional parameter `weights` is used to determine the weights of the data points. It
can be either a single number or a one-dimentional array of length `N`. If a number is
given, that value is used as the error (in magnitude units) for each data point. If an
array is given, it is assumed to contain the error for each data point. If no value is
given, a default error of 0.03 mag is used for each data point.

The output value `fit` is an array containing the best-fit values of H, G1 and G2.

The output value `covmatrix` is the covariance matrix of the least squares fit, and can
be used to produce an error estimate (see below).

### Fitting H, G12

The function `fit_HG12` performs the H, G12 fit. The default form is

`fit, covmatrix = fit_HG12(data, weights=None, degrees=True)`

where the parameters work exactly as for H, G1, G2.


### Error estimation

The function `simulate_errors` returns a Monte Carlo error estimate for the parameters,
both for the H,G1,G2 and the H,G12 systems.

`error_estimates = simulate_errors(fit, covmatrix, N=100000)`

The parameters `fit` and `covmatrix` are the best-fit parameter vector and covariance
matrix returned by the fit function.

The optional parameter `N` gives the number of samples for the Monte Carlo estimate. The
default if 100 000.

The output is an array of eight lines and 3 or 2 columns (depending on the number of
parameters). The first line gives the mean of the Monte Carlo sample, the second line
gives the median. The six lines after that give first the lower and then the upper
percentiles for one, two and three sigma confidence.
