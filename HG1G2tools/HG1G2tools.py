#!/usr/bin/env python
#coding: utf8


from .piecewisefunctions import *
from .spline import Spline
from .HG1G2conversions import *
from numpy import array, radians, matrix, zeros
from numpy.linalg import lstsq as solve_leastsq

class Basis:
    """Container of basis functions for HG1G2 system."""
    def __init__(self, filename=None):
        """Initialize the basis.
        
        The filename parameter is for future purposes, and using it
        will cause a NotImplementedError. Without parameters, the default
        basis is initialized."""
        if filename:
            raise NotImplementedError("Only the default basis is implemented so far.")
        # Construct default basis
        # Define constants for brevity
        r7 = radians(7.5)
        r30 = radians(30)
        r150 = radians(150.0)
        # Make the basis functions
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
             [0.925,0.628842,0.317555,0.127164,0.0223739,0.000165057],
             [-0.572958, -8.6573138e-8]
        )
        # Third basis function, spline and constant zero
        phi3[0, r30] = Spline(
            [radians(x) for x in [0.0, 0.3, 1.0, 2.0, 4.0, 8.0, 12.0, 20.0, 30.0]],
            [1.,0.833812,0.577354,0.421448,0.231742,0.103482,0.0617335,0.016107,0.0],
            [-0.106301, 0.0]
        )
        phi3[r30, r150] = 0.0
        self.version = "20101000"
    def __call__(self, i, x):
        """Return value of the ith basis function at x."""
        return self.functions[i](x)
    def fit_HG1G2(self, data, weight=None, degrees=True):
        """Fit HG1G2 system to data using this basis."""
        Ndata = data.shape[0]
        Ncols = data.shape[1]
        Nfuncs = len(self.functions)
        # Convert angles to radians if needed
        if degrees:
            xval = radians(data[:,0])
        else:
            xval = data[:,0]
        yval = 10**(-0.4 * data[:,1])
        if weight:
            if isinstance(weight, Number):
                errors = zeros(Ndata) + weight
            else:
                errors = weight
        else:
            errors = zeros(Ndata)+0.03
        sigmas = yval * (10**(0.4*errors) - 1)
        Amatrix = matrix(zeros((Ndata, Nfuncs)))
        for i in xrange(Ndata):
            for j in xrange(Nfuncs):
                Amatrix[i,j] = self(j, xval[i]) / sigmas[i]
        yval = 1/(10**(0.4*errors) - 1)
        params = solve_leastsq(Amatrix, yval)
        return a1a2a3_to_HG1G2(params[0])
        

def form_base(filename=None):
    """Return a basis object."""
    return Basis(filename)

def fit_HG1G2(basis):
    pass
