import numpy as np

from .fitzpatrick99 import Fitzpatrick99
from .gordon import Gordon03_SMCBar
from ..external.ezunits import unit


def testunit():
    # check that things look correct
    # -> make some plots
    import pylab

    x = np.arange(0.1, 10, 0.1)   # in um^-1
    lamb = 1.e4 / x * unit['angstrom']

    # ccm  = Cardelli()
    f99  = Fitzpatrick99()
    gsmc = Gordon03_SMCBar()

    fig = pylab.figure()
    plt = fig.add_subplot(1, 1, 1)

    Rv_vals = np.arange(2, 6, dtype=float)
    for Rv in Rv_vals:
        # yccm = ccm.function(lamb, Rv=Rv)
        yf99 = f99.function(lamb, Rv=Rv)

        # pylab.plot(x,yccm,label='CCM, Rv=%0.1f' % (Rv) )
        plt.plot(x, yf99, label='F99, Rv=%0.1f' % (Rv) )

    ygsmc = gsmc.function(lamb)
    plt.plot(x, ygsmc, label='G. SMC')

    mixlaw = f99 + gsmc
    ymix = mixlaw(lamb, Rv=3.1, f_bump=0.75)
    plt.plot(x, ymix, label='Mixture f(bump)=0.75')

    ymix = mixlaw(lamb, Rv=3.1, f_bump=0.5)
    plt.plot(x, ymix, label='Mixture f(bump)=0.5')

    ymix = mixlaw(lamb, Rv=3.1, f_bump=0.25)
    plt.plot(x, ymix, label='Mixture f(bump=0.25')

    plt.set_ylabel('A($\lambda$)/A(V)')
    plt.set_xlabel('1/x [$\mu$m$^{-1}$]')

    plt.legend(loc=0, frameon=False)

    pylab.show()

if __name__ == "__main__":
    testunit()
