HG1G2 tools
===========

This is a Python package for fitting the H,G1,G2 magnitude system to
asteroid photometric data.

This package is still under heavy developement and anything is
susceptible to change until further notice.

Installation
------------

Running ``python setup.py install`` should be enough. NumPy is required.

Usage
-----

::

    import HG1G2tools

    data = loadtxt("44_Nysa.dat", skiprows=1)
    fit, covmatrix = HG1G2tools.fit_HG1G2(data)
    error_estimates = HG1G2tools.simulate_errors()

