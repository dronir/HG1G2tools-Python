    
from __future__ import division

class Spline:
    def __init__(self, xval, yval, deriv):
        N = len(xval)
        A = [0.0 for i in xrange(N)]
        B = [1.0 for i in xrange(N)]
        C = [0.0 for i in xrange(N)]
        R = [0.0 for i in xrange(N)]
        gamma = [0.0 for i in xrange(N)]
        U = [0.0 for i in xrange(N)]
        R[0] = deriv[0]
        for i in xrange(1,N-1):
            A[i] = 1 / (xval[i] - xval[i-1])
            B[i] = 2/(xval[i] - xval[i-1]) + 2/(xval[i+1]-xval[i])
            C[i] = 1/(xval[i+1] - xval[i])
            R[i] += 3 * (yval[i]-yval[i-1])/(xval[i]-xval[i-1])**2
            R[i] += 3 * (yval[i+1]-yval[i])/(xval[i+1]-xval[i])**2
        R[N-1] = deriv[1]
        
        beta = B[0]
        U[0] = R[0]/beta
        for i in xrange(1,N):
            gamma[i] = C[i-1] / beta
            beta = B[i] - A[i] * gamma[i]
            U[i] = (R[i] - A[i]*U[i-1]) / beta
        
        for i in xrange(N-2,-1,-1):
            U[i] = U[i] - gamma[i+1]*U[i+1]
        self.knots = xval[:]
        self.yval = yval[:]
        self.deriv = U[:]
    def __call__(self, x):
        i = None
        for j in xrange(len(self.knots)-1):
            if self.knots[j] <= x and self.knots[j+1] > x:
                i = j
                break
        if i==None:
            raise ValueError("Value %f not in spline support." % x)
        x1 = self.knots[i]
        x2 = self.knots[i+1]
        y1 = self.yval[i]
        y2 = self.yval[i+1]
        d1 = self.deriv[i]
        d2 = self.deriv[i+1]
        xd = x2-x1
        yd = y2-y1
        a = d1*xd - yd
        b = -d2*xd + yd
        t = (x-x1)/xd
        y = (1-t)*y1 + t*y2 + t*(1-t)*(a*(1-t)+b*t)
        return y
