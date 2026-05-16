# Cole-Cole plot for SLS model 2
import numpy as np
import matplotlib.pyplot as plt

k = 10
param_text = '(k = {0})'.format(k)

xi_min, xi_max = -2, 2

xi = np.logspace(xi_min, xi_max, 200)
strMod = [1/k + (1 - 1/k)*x**2 / (1 + x**2) for x in xi]
losMod = [(1 - 1/k)*x / (1 + x**2) for x in xi]

fig = plt.figure(figsize=(8,8), tight_layout=True)
ax = fig.add_subplot(111)
ax.set_title('Cole-Cole plot for SLS model II '+param_text)
ax.set_xlim(-0.05*np.max(strMod), 1.05*np.max(strMod))
ax.set_ylim(-0.05*np.max(strMod), 1.05*np.max(strMod))
ax.set_xlabel("E' / Eins")
ax.set_ylabel('E" / Eins')
ax.scatter(strMod, losMod)
ax.grid()
ax.set_axisbelow(True)

fig.savefig('./png/SLS2_cole_cole.png')

plt.show()