import numpy as np
from scipy import interpolate

from .extinction import ExtinctionLaw, val_in_unit


class Fitzpatrick99(ExtinctionLaw):
    """
    Fitzpatrick (1999, PASP, 111, 63) [1999PASP..111...63F]_
    R(V) dependent extinction curve that explicitly deals with optical/NIR
    extinction being measured from broad/medium band photometry.
    Based on fm_unred.pro from the IDL astronomy library


    .. [1999PASP..111...63F] http://adsabs.harvard.edu/abs/1999PASP..111...63F
    """
    def __init__(self):
        self.name = 'Fitzpatrick99'

    def function(self, lamb, Av=1, Rv=3.1, Alambda=True, **kwargs):
        """
        Fitzpatrick99 extinction Law
        Lamb is input in Anstroms

        Parameters
        ----------
        lamb: float or ndarray(dtype=float)
            wavelength [in Angstroms] at which evaluate the law.

        Av: float
            desired A(V) (default 1.0)

        Rv: float
            desired R(V) (default 3.1)

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

        c2 = -0.824 + 4.717 / Rv
        c1 = 2.030 - 3.007 * c2
        c3 = 3.23
        c4 = 0.41
        x0 = 4.596
        gamma = 0.99

        x = 1.e4 / _lamb
        k = np.zeros(np.size(x))

        # compute the UV portion of A(lambda)/E(B-V)
        xcutuv = 10000.0 / 2700.
        xspluv = 10000.0 / np.array([2700., 2600.])
        yspluv = c1 + (c2 * xspluv) + c3 * ((xspluv) ** 2) / ( ((xspluv) ** 2 - (x0 ** 2)) ** 2 + (gamma ** 2) * ((xspluv) ** 2 ))
        ind = (x >= xcutuv)

        if True in ind:
            k[ind] = c1 + (c2 * x[ind]) + c3 * ((x[ind]) ** 2) / ( ((x[ind]) ** 2 - (x0 ** 2)) ** 2 + (gamma ** 2) * ((x[ind]) ** 2 ))

            # FUV portion
            fuvind = np.where(x >= 5.9)
            k[fuvind] += c4 * (0.5392 * ((x[fuvind] - 5.9) ** 2) + 0.05644 * ((x[fuvind] - 5.9) ** 3))

            k[ind] += Rv
            yspluv += Rv

        # Optical/NIR portion

        ind = x < xcutuv
        if True in ind:
            xsplopir = np.zeros(7)
            xsplopir[0] = 0.0
            xsplopir[1: 7] = 10000.0 / np.array([26500.0, 12200.0, 6000.0, 5470.0, 4670.0, 4110.0])

            ysplopir = np.zeros(7)
            ysplopir[0: 3] = np.array([0.0, 0.26469, 0.82925]) * Rv / 3.1

            ysplopir[3: 7] = np.array([np.poly1d([2.13572e-04, 1.00270, -4.22809e-01])(Rv),
                                       np.poly1d([-7.35778e-05, 1.00216, -5.13540e-02])(Rv),
                                       np.poly1d([-3.32598e-05, 1.00184,  7.00127e-01])(Rv),
                                       np.poly1d([ 1.19456, 1.01707, -5.46959e-03, 7.97809e-04, -4.45636e-05][::-1])(Rv)])

            tck = interpolate.splrep(np.hstack([xsplopir, xspluv]), np.hstack([ysplopir, yspluv]), k=3)
            k[ind] = interpolate.splev(x[ind], tck)

        # convert from A(lambda)/E(B-V) to A(lambda)/A(V)
        k /= Rv

        if (Alambda):
            return(k * Av)
        else:
            return(k * Av * (np.log(10.) * 0.4))
