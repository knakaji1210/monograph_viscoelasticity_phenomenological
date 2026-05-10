# Cole-Cole plot for Maxwell model
import numpy as np
import matplotlib.pyplot as plt

xi_min, xi_max = -2, 2

xi = np.logspace(xi_min, xi_max, 200)
strMod = np.array([x**2 / (1 + x**2) for x in xi])
losMod = np.array([x / (1 + x**2) for x in xi])

fig = plt.figure(figsize=(8,8), tight_layout=True)
ax = fig.add_subplot(111)
ax.set_title('Cole-Cole plot for Maxwell model')
ax.set_xlim(-0.05*np.max(strMod), 1.05*np.max(strMod))
ax.set_ylim(-0.05*np.max(strMod), 1.05*np.max(strMod))
ax.set_xlabel(r'$E^{\prime}$ / $E_i$')
ax.set_ylabel(r'$E^{{\prime\prime}}$ / $E_i$')
ax.scatter(strMod, losMod)
ax.grid()
ax.set_axisbelow(True)

savefile = './png/Maxwell_cole_cole.png'
fig.savefig(savefile, dpi=300)

plt.show()