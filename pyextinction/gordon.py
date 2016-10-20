import numpy as np
from scipy import interpolate

from .extinction import ExtinctionLaw, val_in_unit


class Gordon03_SMCBar(ExtinctionLaw):
    """ Gordon et al. 2003 (ApJ, 594:279-293)

    In principle, :math:`R_V` has no impact on this law: according to Gordon et al
    (2003), the average value of :math:`R_V` is fixed to :math:`2.74 \pm 0.13`
    In practice we offer to change this value, if one wants to explore the
    uncertainties.

    Attributes
    ----------
    Rv: float
        desired default R(V), can be replaced during calling sequences
    """
    def __init__(self, Rv=2.74):
        """
        Parameters
        ----------
        Rv: float
            desired R(V) (default internal value given at initialization)
        """
        self.name = 'Gordon et al. 2003 SMCBar'
        self.Rv = Rv

    def function(self, lamb, Av=1, Rv=None, Alambda=True,  **kwargs):
        """
        Lamb is input in Anstroms
        Note that Rv is not given as a variable in the paper of reference

        Parameters
        ----------
        lamb: float or ndarray(dtype=float)
            wavelength [in Angstroms] at which evaluate the law.

        Av: float
            desired A(V) (default 1.0)

        Rv: float
            desired R(V) (default internal value given at initialization)

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

        if Rv is None:
            Rv = self.Rv

        c1 = -4.959 / Rv
        c2 = 2.264 / Rv
        c3 = 0.389 / Rv
        c4 = 0.461 / Rv
        x0 = 4.6
        gamma = 1.0

        x = 1.e4 / _lamb
        k = np.zeros(np.size(x))

        # UV part
        xcutuv = 10000.0 / 2700.
        xspluv = 10000.0 / np.array([2700., 2600.])
        yspluv = 1.0 + c1 + (c2 * xspluv) + c3 * ((xspluv) ** 2) / ( ((xspluv) ** 2 - (x0 ** 2)) ** 2 + (gamma ** 2) * ((xspluv) ** 2 ))

        ind = np.where(x >= xcutuv)
        if np.size(ind) > 0:
            k[ind] = 1.0 + c1 + (c2 * x[ind]) + c3 * ((x[ind]) ** 2) / ( ((x[ind]) ** 2 - (x0 ** 2)) ** 2 + (gamma ** 2) * ((x[ind]) ** 2 ))

            ind = np.where(x >= 5.9)
            k[ind] += c4 * (0.5392 * ((x[ind] - 5.9) ** 2) + 0.05644 * ((x[ind] - 5.9) ** 3))

        # Opt/NIR part
        ind = np.where(x < xcutuv)
        if np.size(ind) > 0:
            xsplopir = np.zeros(9)
            xsplopir[0] = 0.0
            xsplopir[1: 10] = 1.0 / np.array([2.198, 1.65, 1.25, 0.81, 0.65, 0.55, 0.44, 0.37])

            # Values directly from Gordon et al. (2003)
            # ysplopir =  np.array([0.0,0.016,0.169,0.131,0.567,0.801,1.00,1.374,1.672])
            # K & J values adjusted to provide a smooth, non-negative cubic spline interpolation
            ysplopir = np.array([0.0, 0.11, 0.169, 0.25, 0.567, 0.801, 1.00, 1.374, 1.672])

            tck = interpolate.splrep(np.hstack([xsplopir, xspluv]), np.hstack([ysplopir, yspluv]), k=3)
            k[ind] = interpolate.splev(x[ind], tck)

        if (Alambda):
            return(k * Av)
        else:
            return(k * Av * (np.log(10.) * 0.4 ))
