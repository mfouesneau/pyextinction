"""
Dust Extinction laws
--------------------


.. note::

    This module is able to handle values with units
"""
from .helpers import val_in_unit, isNestedInstance

__version__ = '1.0'
__all__ = ['ExtinctionLaw', 'MixtureLaw']


class ExtinctionLaw(object):
    """ Template class """
    def __init__(self):
        self.name = 'None'

    def __repr__(self):
        return '{0:s}\n{1:s}'.format(self.name, object.__repr__(self))

    def function(self, lamb, *arg, **kwargs):
        """ expected to contain a function of lambda that return the
        extinction values

        Parameters
        ----------
        lamb: ndarray
            wavelength

        Returns
        -------
        val: ndarray
            expected values of the law evaluated at lamb
        """
        raise NotImplementedError

    def __call__(self, *args, **kwargs):
        """ Make the extinction law callable object using :func:`self.function`"""
        return self.function(*args, **kwargs)

    def isvalid(self, *args, **kwargs):
        """ Check if the current arguments are in the validity domain of the law
        Must be redefined if any restriction applies to the law
        """
        return True

    def __add__(self, other):
        return MixtureLaw(A=self, B=other)


class MixtureLaw(ExtinctionLaw):
    """
    Mixture of extinction laws allowing to vary for instance the bump amplitude
    in the extinction law

    ..math::

            f_A * A(*args, **kwargs) + (1 - f_A) * B(*args, **kwargs)

    .. example::

        l = Fitzpatrick99() + Gordon03_SMCBar()

    """
    def __init__(self, A=None, B=None, name=None):
        """ Constructor

        Parameters
        ----------
        A: ExtinctionLaw
            Component A

        B: ExtinctionLaw
            Component B
        """
        if not (isNestedInstance(A, ExtinctionLaw) &
                isNestedInstance(B, ExtinctionLaw)):
            raise ValueError('Expecting ExtinctionLaw instances')
        self.A = A
        self.B = B
        self.name = name or '(' + self.A.name + ', ' + self.B.name + ')'

    def function(self, lamb, Av=1, Rv_A=None, Alambda=True, f_A=0.5, Rv_B=None,
                 Rv=None, **kwargs):
        """
        Lamb as to be in Angstroms!!!

        Parameters
        ----------
        lamb: float or ndarray(dtype=float)
            wavelength [in Angstroms] at which evaluate the law.

        Av: float
            desired A(V) (default 1.0)

        Alambda: bool
            if set returns +2.5*1./log(10.)*tau, tau otherwise

        f_A: float
            set the mixture ratio between the two laws (default 0.5)

        Rv_A: float
            extinction param. on the Law A

        Rv_B: float
            extinction param. on the bumpless component

        Rv: float
            effective R(V) according to the mixture

        Returns
        -------
        r: float or ndarray(dtype=float)
            attenuation as a function of wavelength
            depending on Alambda option +2.5*1./log(10.)*tau,  or tau

        .. math::

            f_A * A(*args, **kwargs) + (1. - f_A) * B(*args, **kwargs)
        """
        u_lamb = val_in_unit('lamb', lamb, 'angstrom')

        if Rv_A is None:
            Rv_A = getattr(self.A, 'Rv', None)

        if Rv_B is None:
            Rv_B = getattr(self.B, 'Rv', None)

        if sum([Rv_A is None, Rv_B is None, Rv is None]) >= 2:
            raise ValueError('Must provide at least 2 Rv values')

        if Rv_A is None:
            Rv_A = self.get_Rv_A(Rv, f_A, Rv_B)
        if Rv_B is None:
            Rv_B = self.get_Rv_B(Rv, Rv_A, f_A)

        return (f_A * self.A.function(u_lamb, Av=Av, Rv=Rv_A, Alambda=Alambda)
                + (1. - f_A) * self.B.function(u_lamb, Av=Av, Alambda=Alambda,
                                               Rv=Rv_B)
                )

    def isvalid(self, Av=None, Rv=None, f_A=0.5, Rv_A=None, Rv_B=None):
        """ Test the validity of an extinction vector (Av, Rv, Rv_A, Rv_B, fbump)

        .. math::
            Law = f_A * A(lamb, Av=Av, Rv=Rv_A) + (1. - f_A) * B(lamb, Av=Av, Rv=Rv_B)

        The validity impose :math:`R_V` ranges and  to be a fraction (i.e.,
        between 0 and 1)

        At least 2 out of the 3 :math:`R_V` values must be provided, the 3rd
        will be computed if missing.

        Parameters
        ----------
        Av: float
            Av value (any value is allowed, even <0)

        Rv, Rv_A, Rv_B: float, float, float
            effective Rv, A component and B component Rv values, respectively.
            At least 2 must be provided.

        f_A: float
            Mixture ratio between the two components

        Returns
        -------
        r: bool
            True, if the values a coherent with the definition.
        """

        if Rv_B is None and hasattr(self.B, 'Rv'):
            Rv_B = self.B.Rv

        if Rv_A is None and hasattr(self.A, 'Rv'):
            Rv_A = self.A.Rv

        # if we do not have at least 2 of the 3 Rvs defined then it's invalid
        if sum([Rv_A is None, Rv_B is None, Rv is None]) >= 2:
            return False

        if Rv_A is None:
            Rv_A = self.get_Rv_A(Rv, f_A, Rv_B=Rv_B)
        if Rv is None:
            Rv = self.get_Rv(Rv_A, f_A, Rv_B=Rv_B)
        if Rv_B is None:
            Rv_B = self.get_Rv_B(Rv, Rv_A, f_A)

        # f_A is a fraction and any Rv is limited to [2.0, 6.0]
        return ((0. <= f_A <= 1.) & (2.0 <= Rv_B <= 6.0) &
                (2.0 <= Rv_A <= 6.0) & (2.0 <= Rv <= 6.0))

    def get_Rv_A(self, Rv, f_A=0.5, Rv_B=None):
        """ Returns the equivalent Rv to use in the bump component
            Law = f_A * A (lamb, Av=Av, Rv=Rv_A) + (1. - f_A) * B(lamb, Av=Av, Rv=Rv_B)

            and Rv_A is such that:

            ..math::

                1 / Rv = f_A / Rv_A + (1 - f_A) / Rv_B

                Rv_A = 1. / (1. / (Rv * f_A) - (1. - f_A) / (f_A * Rv_B))

            not that Gordon03_SMCBar has a fixed Rv=2.74
        """
        if Rv_B is None and hasattr(self.B, 'Rv'):
            Rv_B = self.B.Rv

        return 1. / (1. / (Rv * f_A) - (1. - f_A) / (f_A * Rv_B))

    def get_Rv(self, Rv_A=None, f_A=0.5, Rv_B=None):
        """ Returns the equivalent effective Rv according to the mixture

        ..math::

            Law = f_A * A (lamb, Av=Av, Rv=Rv_A) + (1. - f_A) * B(lamb, Av=Av, Rv=Rv_B)

        and Rv is such that:

        ..math::

            1 / Rv = f_A / Rv_A + (1 - f_A) / Rv_B

            Rv_A = 1. / (1. / (Rv * f_A) - (1. - f_A) / (f_A * Rv_B))
        """
        if Rv_B is None and hasattr(self.B, 'Rv'):
            Rv_B = self.B.Rv

        if Rv_A is None and hasattr(self.A, 'Rv'):
            Rv_A = self.A.Rv

        return 1. / (f_A / Rv_A + (1 - f_A) / Rv_B)

    def get_Rv_B(self, Rv, Rv_A=None, f_A=0.5):
        """ Returns the equivalent Rv to use in the bumpless component

        .. math::
            Law = f_A * A (lamb, Av=Av, Rv=Rv_A) + (1. - f_A) * B(lamb, Av=Av, Rv=Rv_B)

        and Rv_B is such that

        .. math::

            1 / Rv = f_A / Rv_A + (1 - f_A) / Rv_B

            Rv_A = 1. / (1. / (Rv * f_A) - (1. - f_A) / (f_A * Rv_B))
        """
        if Rv_A is None and hasattr(self.A, 'Rv'):
            Rv_A = self.A.Rv

        return (1. - f_A) / (1. / Rv - f_A / Rv_A)
