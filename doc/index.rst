.. pyextinction documentation master file, created by
   sphinx-quickstart on Thu Oct 20 10:29:20 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pyextinction - A tool using stellar extinction curves
=====================================================

In many applications, one need physical models that predict spectra or SEDs of
star extinguished by dust.  
Interstellar dust extinguishes stellar light as it
travels from the star’s surface to the observer. The wavelength dependence of
the extinction from the UV to the NIR has been measured along many sightlines in
the Milky Way (Cardelli et al. 1989; Fitzpatrick 1999; Valencic et al. 2004;
Gordon et al. 2009) and for a handful of sightlines in the Magellanic Clouds
(Gordon & Clayton 1998; Misselt et al. 1999; Maız Apellaniz & Rubio 2012) as
well as in M31 (Bianchi et al. 1996, Clayton et al. 2015, submitted).

The observations show a wide range of dust column normalized extinction curves,
:math:`A(\lambda) / A(V)`.  This package provides a common interface to many
commonly used extinction laws

Package main content
~~~~~~~~~~~~~~~~~~~~

* :class:`pyextinction.Cardelli`, Cardelli, Clayton, and Mathis (1989, ApJ, 345, 245)
* :class:`pyextinction.Calzetti`, Calzetti et al. (2000, ApJ 533, 682)
* :class:`pyextinction.Fitzpatrick`, Fitzpatrick (1999, PASP, 111, 63) 
* :class:`pyextinction.Gordon03_SMCBar`, Gordon et al. 2003 (ApJ, 594:279-293)

Once could also combine laws into a single one. (see below)

Contents:

.. toctree::
   :maxdepth: 2

   modules

Quick Start
~~~~~~~~~~~

.. code-block:: python 

        import pylab as plt
        import pyextinction
        from pyextinction import unit
        import numpy as np


        # list of law to test
        laws = (pyextinction.Cardelli(),
                pyextinction.Fitzpatrick99(),
                pyextinction.Gordon03_SMCBar()
                )

        #define the wave numbers
        x = np.arange(0.1, 10, 0.1)   # in um^-1
        lamb = (1e4 / x) * unit['angstrom']

        Rv = 3.1
        for l in laws:
                plt.plot(x, l(lamb, Rv=Rv), label=l.name, lw=2)

        plt.legend(loc='upper left', frameon=False)
        plt.xlabel(r'Wave number [$\mu$m$^{-1}$]')
        plt.ylabel(r'$A(\lambda) / A(V)$')


.. image:: multiple_laws.png
        :scale: 70 %


.. code-block:: python


        Rv_vals = np.arange(2, 6, dtype=float)
        l = pyextinction.Fitzpatrick99()

        for Rv in Rv_vals:
                plt.plot(x, l(lamb, Rv=Rv), label='Rv={0:0.1f}'.format(Rv), lw=2)
        plt.legend(loc='upper left', frameon=False)
        plt.xlabel(r'Wave number [$\mu$m$^{-1}$]')
        plt.ylabel(r'$A(\lambda) / A(V)$')


.. image:: multiple_Rv.png
        :scale: 70 %

Mixing different laws
~~~~~~~~~~~~~~~~~~~~~

we introduce a mixture model with two components $A$ and $B$ to describe the
full range of observed extinction curves in the Local Group.

When $f_A$ gives the fraction of the :math:`A`-type extinction and :math:`(1 − f_A)` the
fraction of the :math:`B`-type extinction, this mixture is defined as follow: 

.. math::

        \frac{A(\lambda)}{A(V)} = f_A\,\left[\frac{A(\lambda)}{A(V)}\right]_A + (1 - f_A)\,\left[\frac{A(\lambda)}{A(V)}\right]_B,

The effective :math:`R(V)` of the mixture extinction curve model is then given by:

.. math::

        R(V)^{-1} = f_A\,R_A(V)^{-1} + (1 - f_A)\,R_B(V)^{-1}.


This code is made flexible enough that any laws can be combined. 

Below we use `Fitzpatrick` (A) definition with the `SMCBar` (B) with fix
:math:`R_B(V) = 2.74` (Gordon et al. 2003), as we did in Gordon et al.
2016. 
The range of observed :math:`R_A(V)` is between :math:`2.0` and :math:`6.0` and
this results in the parameter space defined by :math:`(R(V), f_A)` not being
completely filled.

.. code-block:: python

        import pyextinction
        mixture = pyextinction.Fitzpatrick99() + pyextinction.Gordon03_SMCBar()

        Rv = 3.1
        f_A_vals = (0.1, 0.25, 0.5, 0.75, 0.9)

        for f_A in f_A_vals:
                plt.plot(x, mixture(lamb, Rv=Rv, f_A=f_A),
                         label=r'f$_A$={0:0.2f}'.format(f_A), lw=2)
        plt.legend(loc='upper left', frameon=False, bbox_to_anchor=(1.05, 1.05) )
        plt.xlabel(r'Wave number [$\mu$m$^{-1}$]')
        plt.ylabel(r'$A(\lambda) / A(V)$')
        figrc.hide_axis('top right'.split())

.. image:: mixture.png
        :scale: 70 %



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

