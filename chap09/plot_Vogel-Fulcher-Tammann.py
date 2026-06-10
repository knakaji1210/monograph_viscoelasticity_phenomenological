# Plot of Vogel-Fulcher-Tammann equation
# 260322修正（講義：高分子レオロジーのため）

import numpy as np
import matplotlib.pyplot as plt


temp_glass = 200
temp_vogel = temp_glass - 50

temp_A = 1727

preF = 34.5

eta_inf = 10**(-3) # unit Pa s

temp_min, temp_max = temp_vogel, 300

temp = np.linspace(temp_min, temp_max, 50)
Arrhenius = [eta_inf*np.exp(preF*temp_glass/T) for T in temp]
VFT = [eta_inf*np.exp(temp_A/(T - temp_vogel)) for T in temp]
log_Arrhenius = [np.log10(arr) for arr in Arrhenius]
log_VFT = [np.log10(vft) for vft in VFT]
func_min = np.min(log_Arrhenius)
func_max = np.max(log_Arrhenius)
ylim = [0.1*func_min, func_max*2]

fig = plt.figure(figsize=(8,5), tight_layout=True)
ax = fig.add_subplot(111)
ax.set_title('Arrhnius vs Vogel-Fulcher-Tammann')
ax.set_xlim(temp_min-10, temp_max+10)
ax.set_ylim(ylim[0], ylim[1])
ax.set_xlabel(r'$T$ /K')
ax.set_ylabel(r'log($\eta$($T$) /Pa s)')
ax.scatter(temp, log_Arrhenius, label='Arrhenius')
ax.scatter(temp, log_VFT, label='Vogel-Fulcher-Tammann')

ax.vlines([temp_glass], ylim[0], ylim[1], 'g', ls='--')
ax.vlines([temp_vogel], ylim[0], ylim[1], 'b', ls='--')
ax.hlines([12], temp_min-10, temp_glass, 'r', ls='--')

ax.legend(loc='upper right') 
ax.grid()
ax.set_axisbelow(True)

fig_text1 = r'$T_{{V}}$ = {0:.0f} K'.format(temp_vogel) 
fig_text2 = r'$T_{{g}}$ = {0:.0f} K'.format(temp_glass) 
fig.text(0.18, 0.15, fig_text1, c='b')
fig.text(0.41, 0.15, fig_text2, c='g')
fig.text(0.09, 0.36, "12", c='r')

fig.savefig('./png/VFT.png', dpi=300)

plt.show()