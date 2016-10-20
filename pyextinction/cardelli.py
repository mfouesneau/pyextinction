import numpy as np
from .extinction import ExtinctionLaw, val_in_unit


class Cardelli(ExtinctionLaw):
    """ Cardelli, Clayton, and Mathis (1989, ApJ, 345, 245)"""
    def __init__(self):
        self.name = 'Cardelli'

    def function(self, lamb, Av=1., Rv=3.1, Alambda=True, **kwargs):
        """
        Cardelli extinction Law

        Parameters
        ----------
        lamb: float or ndarray(dtype=float)
            wavelength [in Angstroms] at which evaluate the law.

        Av: float
            desired A(V) (default: 1.0)

        Rv: float
            desired R(V) (default: 3.1)

        Alambda: bool
            if set returns +2.5*1./log(10.)*tau, tau otherwise

        Returns
        -------
        r: float or ndarray(dtype=float)
            attenuation as a function of wavelength
            depending on Alambda option +2.5*1./log(10.)*tau,  or tau
        """
        _lamb = val_in_unit('lamb', lamb, 'angstrom').magnitude

        if isinstance(_lamb, float) or isinstance(_lamb, np.float_):
            _lamb = np.asarray([_lamb])
        else:
            _lamb = _lamb[:]

        # init variables
        x = 1.e4 / _lamb  # wavenumber in um^-1
        a = np.zeros(np.size(x))
        b = np.zeros(np.size(x))
        # Infrared (Eq 2a,2b)
        ind = np.where((x >= 0.3) & (x < 1.1))
        a[ind] =  0.574 * x[ind] ** 1.61
        b[ind] = -0.527 * x[ind] ** 1.61
        # Optical & Near IR
        # Eq 3a, 3b
        ind = np.where((x >= 1.1) & (x <= 3.3))
        y = x[ind] - 1.82
        a[ind] = 1. + 0.17699 * y - 0.50447 * y ** 2 - 0.02427 * y ** 3 + 0.72085 * y ** 4 + 0.01979 * y ** 5 - 0.77530 * y ** 6 + 0.32999 * y ** 7
        b[ind] =      1.41338 * y + 2.28305 * y ** 2 + 1.07233 * y ** 3 - 5.38434 * y ** 4 - 0.62251 * y ** 5 + 5.30260 * y ** 6 - 2.09002 * y ** 7
        # UV
        # Eq 4a, 4b
        ind = np.where((x >= 3.3) & (x <= 8.0))
        a[ind] =  1.752 - 0.316 * x[ind] - 0.104 / ((x[ind] - 4.67) ** 2 + 0.341)
        b[ind] = -3.090 + 1.825 * x[ind] + 1.206 / ((x[ind] - 4.62) ** 2 + 0.263)

        ind = np.where((x >= 5.9) & (x <= 8.0))
        Fa     = -0.04473 * (x[ind] - 5.9) ** 2 - 0.009779 * (x[ind] - 5.9) ** 3
        Fb     =  0.21300 * (x[ind] - 5.9) ** 2 + 0.120700 * (x[ind] - 5.9) ** 3
        a[ind] = a[ind] + Fa
        b[ind] = b[ind] + Fb
        # Far UV
        # Eq 5a, 5b
        ind = np.where((x >= 8.0) & (x <= 10.0))
        # Fa = Fb = 0
        a[ind] = -1.073 - 0.628 * (x[ind] - 8.) + 0.137 * ((x[ind] - 8.) ** 2) - 0.070 * (x[ind] - 8.) ** 3
        b[ind] = 13.670 + 4.257 * (x[ind] - 8.) + 0.420 * ((x[ind] - 8.) ** 2) + 0.374 * (x[ind] - 8.) ** 3

        # Case of -values x out of range [0.3,10.0]
        ind = np.where((x > 10.0) | (x < 0.3))
        a[ind] = 0.0
        b[ind] = 0.0

        # Return Extinction vector
        # Eq 1
        if (Alambda):
            return( ( a + b / Rv ) * Av)
        else:
            # return( 1./(2.5 * 1. / np.log(10.)) * ( a + b / Rv ) * Av)
            return( 0.4 * np.log(10.) * ( a + b / Rv ) * Av)
