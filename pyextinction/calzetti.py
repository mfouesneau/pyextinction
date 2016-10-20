import numpy as np
from .extinction import ExtinctionLaw, val_in_unit


class Calzetti(ExtinctionLaw):
    """
    Calzetti et al. (2000, ApJ 533, 682) developed a recipe for dereddening the
    spectra of galaxies where massive stars dominate the radiation output,
    strictly valid between 0.12 to 2.2 microns, and extrapolated from 0.12 down
    to 0.0912 microns.

    Note that the supplied color excess should be that derived for the
    stellar  continuum, :math:`E_{B-V}(stars)`, which is related to the reddening
    derived from the gas, :math:`E_{B-V}(gas)`, via the Balmer decrement by
    :math:`E_{B-V}(stars) = 0.44 \\times E_{B-V}(gas)`

    :math:`R_V` - Ratio of total to selective extinction, default is 4.05.
    Calzetti et al. (2000) estimate :math:`R_V = 4.05 \pm 0.80` from optical-IR
    observations of 4 starbursts.
    """
    def __init__(self):
        self.name = 'Calzetti'

    def function(self, lamb, Av=1, Rv=4.05, Alambda=True, **kwargs):
        """
        Returns Alambda or tau for a Calzetti law Lamb is input in Angstroms

        Parameters
        ----------
        lamb: float or ndarray(dtype=float)
            wavelength [in Angstroms] at which evaluate the law.

        Av: float
            desired A(V) (default 1.0)

        Rv: float
            desired R(V) (default 4.05)

        Alambda: bool
            if set returns +2.5 * 1. / log(10.) * tau, tau otherwise

        Returns
        -------
        r: float or ndarray(dtype=float)
            attenuation as a function of wavelength
            depending on Alambda option +2.5*1./log(10.)*tau,  or tau
        """
        # handle units
        _lamb = val_in_unit('lamb', lamb, 'angstrom').magnitude

        if isinstance(_lamb, float) or isinstance(_lamb, np.float_):
            _lamb = np.asarray([_lamb])
        else:
            _lamb = _lamb[:]

        _lamb *= 1e-4

        x = 1. / _lamb  # wavenumber in um^-1
        k = np.zeros(np.size(x))

        ind = (_lamb >= 0.630 ) & (_lamb <= 2.2)
        k[ind] = 2.659 * (-1.857 + 1.040 * x[ind]) + Rv

        ind = (_lamb >= 0.0912 ) & (_lamb < 0.630)
        k[ind] = 2.659 * (-2.156 + 1.509 * x[ind] - 0.198 * x[ind] ** 2 + 0.011 * x[ind] ** 3 ) + Rv

        if Alambda:
            return 0.4 * k
        else:
            return 10 ** (0.4 * k)
