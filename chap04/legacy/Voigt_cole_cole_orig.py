# Cole-Cole plot for Kelvin-Voigt model
import numpy as np
import matplotlib.pyplot as plt

xi_min, xi_max = -2, 2

xi = np.logspace(xi_min, xi_max, 200)
strComp = np.array([1 / (1 + x**2) for x in xi])
losComp = np.array([x / (1 + x**2) for x in xi])

fig = plt.figure(figsize=(8,8), tight_layout=True)
ax = fig.add_subplot(111)
ax.set_title('Cole-Cole plot for Voigt model')
ax.set_xlim(-0.05*np.max(strComp), 1.05*np.max(strComp))
ax.set_ylim(-0.05*np.max(strComp), 1.05*np.max(strComp))
ax.set_xlabel(r'$J^{\prime}$ / $J_{{{\infty}}}$')
ax.set_ylabel(r'$J^{{\prime\prime}}$ / $J_{{{\infty}}}$')
ax.scatter(strComp, losComp)
ax.grid()
ax.set_axisbelow(True)

savefile = './png/Voigt_cole_cole.png'
fig.savefig(savefile, dpi=300)

plt.show()