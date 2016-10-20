.. pyextinction documentation master file, created by
   sphinx-quickstart on Thu Oct 20 10:29:20 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pyextinction's documentation!
========================================

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

Once could also combine laws into a single one.

Contents:

.. toctree::
   :maxdepth: 2

   modules

Quick Start
~~~~~~~~~~~


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

