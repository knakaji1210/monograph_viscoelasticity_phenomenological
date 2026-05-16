# Cole-Cole plot for SLS model 1
import numpy as np
import matplotlib.pyplot as plt

k = 10
param_text = '(k = {0})'.format(k)

xi_min, xi_max = -2, 2

xi = np.logspace(xi_min, xi_max, 200)
strComp = [1 + (k - 1) / (1 + x**2) for x in xi]
losComp = [(k - 1)*x / (1 + x**2) for x in xi]

fig = plt.figure(figsize=(8,8), tight_layout=True)
ax = fig.add_subplot(111)
ax.set_title('Cole-Cole plot for SLS model I '+param_text)
ax.set_xlim(-0.05*np.max(strComp), 1.05*np.max(strComp))
ax.set_ylim(-0.05*np.max(strComp), 1.05*np.max(strComp))
ax.set_xlabel("J' / Jinf")
ax.set_ylabel('J" / Jinf')
ax.scatter(strComp, losComp)
ax.grid()
ax.set_axisbelow(True)

fig.savefig('./png/SLS1_cole_cole.png')

plt.show()