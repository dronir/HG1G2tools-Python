#!/usr/bin/env python
#coding: utf8

from .piecewisefunctions import *
from .spline import Spline
from .HG1G2conversions import *
from numpy import array, radians, matrix, zeros, percentile, median
from numpy.linalg import lstsq as solve_leastsq, inv
from numpy.random import multivariate_normal as multinormal

class Basis:
    """Container of basis functions for HG1G2 system."""
    def __init__(self, filename=None):
        """Initialize the basis.
        
        The filename parameter is for future purposes, and using it
        will cause a NotImplementedError. Without parameters, the default
        basis is initialized."""
        if filename:
            raise NotImplementedError("Only the default basis is implemented so far.")

        self.version = "20101000"

        r7 = radians(7.5)
        r30 = radians(30)
        r150 = radians(150.0)

        phi1 = PiecewiseFunction()
        phi2 = PiecewiseFunction()
        phi3 = PiecewiseFunction()
        self.functions = [phi1, phi2, phi3]

        # First basis function, linear and spline
        phi1[0, r7] = lambda x: 1.0 - 1.90985931710274*x
        phi1[r7, r150] = Spline(
            [radians(x) for x in [7.5, 30, 60, 90, 120, 150]],
            [0.75, 0.33486, 0.134106, 0.0511048, 0.0214657, 0.0036397],
            [-1.90986, -0.0913286]
        )

        # Second basis function, linear and spline
        phi2[0, r7] = lambda x: 1.0 - 0.572957795130823*x
        phi2[r7, r150] = Spline(
             [radians(x) for x in [7.5, 30, 60, 90, 120, 150]],
             [0.925, 0.628842, 0.317555, 0.127164, 0.0223739, 0.000165057],
             [-0.572958, -8.6573138e-8]
        )

        # Third basis function, spline and constant zero
        phi3[0, r30] = Spline(
            [radians(x) for x in [0.0, 0.3, 1.0, 2.0, 4.0, 8.0, 12.0, 20.0, 30.0]],
            [1.0, 0.833812, 0.577354, 0.421448, 0.231742, 0.103482, 0.0617335, 0.016107, 0.0],
            [-0.106301, 0.0]
        )
        phi3[r30, r150] = 0.0

    def __len__(self):
        return len(self.functions)
                
    def __call__(self, x):
        """Return values of the basis functions at x."""
        return [self.functions[i](x) for i in xrange(len(self))]
        
    def fit_HG1G2(self, data, weight=None, degrees=True):
        """Fit (H,G1,G2) system to data using this basis."""
        Ndata = data.shape[0]
        Ncols = data.shape[1]

        xval = radians(data[:,0]) if degrees else data[:,0]
        yval = 10**(-0.4 * data[:,1])
        
        if weight:
            if isinstance(weight, Number):
                errors = zeros(Ndata) + weight
            else:
                errors = weight
        else:
            errors = zeros(Ndata)+0.03
        sigmas = yval * (10**(0.4*errors) - 1)
    
        Amatrix = matrix(zeros((Ndata, 3)))
        for i in xrange(Ndata):
            Amatrix[i,:] = array(self(xval[i])) / sigmas[i]
        
        yval = 1/(10**(0.4*errors) - 1)
        params, residual, rank, s = solve_leastsq(Amatrix, yval)
        covMatrix = inv(Amatrix.T * Amatrix)
    
        return a1a2a3_to_HG1G2(params), covMatrix
    
    def fit_HG12(self, data, weight=None, degrees=True):
        """Fit (H,G12) system to data using this basis."""
        Ndata = data.shape[0]
        Ncols = data.shape[1]

        xval = radians(data[:,0]) if degrees else data[:,0]
        yval = 10**(-0.4 * data[:,1])
        
        if weight:
            if isinstance(weight, Number):
                errors = zeros(Ndata) + weight
            else:
                errors = weight
        else:
            errors = zeros(Ndata)+0.03
        sigmas = yval * (10**(0.4*errors) - 1)
        
        pars = (0.7527, 0.06164, -0.9612, 0.6270), (0.9529, 0.02162, -0.6125, 0.5572)
        min_residual = 1e12
        best_params = None
        for parameter_set in pars:
            b1,b0,g1,g0 = parameter_set
            Amatrix = matrix(zeros((Ndata, 2)))
            for i in xrange(Ndata):
                Phi = array(self(xval[i]))
                Gamma1 = b0*Phi[0] + g0*Phi[1] + Phi[2] - b0*Phi[2] - g0*Phi[2]
                Gamma2 = b1*Phi[0] + g1*Phi[1] - (b1+g1)*Phi[2]
                Amatrix[i,0] = Gamma1 / sigmas[i]
                Amatrix[i,1] = Gamma2 / sigmas[i]
            yval = 1/(10**(0.4*errors) - 1)
            params, residual, rank, s = solve_leastsq(Amatrix, yval)
            covMatrix = inv(Amatrix.T * Amatrix)
            if residual[0] < min_residual:
                min_residual = residual[0]
                best_params = params
        return (a1a2_to_HG12(best_params), covMatrix)

def simulate_errors(params, cov, N=100000):
    """Return a parameter error estimate based on the covariance matrix."""
    M = len(params)
    if M == 3:
        convert_from = HG1G2_to_a1a2a3
        convert_to = a1a2a3_to_HG1G2
    elif M == 2:
        convert_from = HG12_to_a1a2
        convert_to = a1a2_to_HG12
    else:
        raise ValueError("Parameter vector length {} (expected 2 or 3).".format(len(params)))
    mean = convert_from(params)
    sample = multinormal(mean, cov, size=N)
    HGsample = array([convert_to(x) for x in sample])
    res = zeros((8, M))
    
    res[0,0:M] = HGsample.mean(axis=0)
    res[1,0:M] = median(HGsample, axis=0)
    res[2,0:M] = percentile(HGsample, 0.13499, axis=0)
    res[3,0:M] = percentile(HGsample, 2.27501, axis=0)
    res[4,0:M] = percentile(HGsample, 15.8655, axis=0)
    res[5,0:M] = percentile(HGsample, 84.1345, axis=0)
    res[6,0:M] = percentile(HGsample, 97.725, axis=0)
    res[7,0:M] = percentile(HGsample, 99.865, axis=0)
    return res

def sample_HG1G2(params, cov, N):
    mean = HG1G2_to_a1a2a3(params)
    sample = multinormal(mean, cov, size=N)
    HGsample = array([a1a2a3_to_HG1G2(x) for x in sample])
    return HGsample

def form_base(filename=None):
    """Return a basis object."""
    return Basis(filename)

def fit_HG1G2(basis):
    pass
