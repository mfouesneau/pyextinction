import figrc
import setup_mpl
setup_mpl.theme()
setup_mpl.solarized_colors()
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
figrc.hide_axis('top right'.split())
plt.tight_layout()
plt.savefig('multiple_laws.png')


plt.figure()

Rv_vals = np.arange(2, 6, dtype=float)
l = pyextinction.Fitzpatrick99()

for Rv in Rv_vals:
    plt.plot(x, l(lamb, Rv=Rv), label='Rv={0:0.1f}'.format(Rv), lw=2)
plt.legend(loc='upper left', frameon=False)
plt.xlabel(r'Wave number [$\mu$m$^{-1}$]')
plt.ylabel(r'$A(\lambda) / A(V)$')
figrc.hide_axis('top right'.split())
plt.tight_layout()
plt.savefig('multiple_Rv.png')


plt.figure()
mixture = pyextinction.Fitzpatrick99() + pyextinction.Gordon03_SMCBar()

Rv = 3.1
f_A_vals = (0.1, 0.25, 0.5, 0.75, 0.9)

for f_A in f_A_vals:
        plt.plot(x, mixture(lamb, Rv=Rv, f_A=f_A),
                 label=r'f$_A$={0:0.2f}'.format(f_A), lw=2)
plt.legend(loc='upper left', frameon=False)
plt.xlabel(r'Wave number [$\mu$m$^{-1}$]')
plt.ylabel(r'$A(\lambda) / A(V)$')
figrc.hide_axis('top right'.split())
plt.tight_layout()
plt.savefig('mixture.png')

