# ordinary differential equation of spring (sinusoidal strain) with animation
'''
バネ要素単独では常微分方程式は不要
'''

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.animation import FuncAnimation

# variables
try:
    E = float(input('modulus [MPa] (default = 0.2 MPa): '))*10**6
except ValueError:
    E = 2*10**5             # [Pa] modulus

l = 0.1                     # [m] equilibrium length
w = 0.5                     # ratio of spring width

# external sinusoidal strain
try:
    eamp = float(input('amplitude for sinusoidal strain [] (default=0.25): '))
except ValueError:
    eamp = 0.25
try:
    T = float(input('peirod for sinusoidal strain [s] (default=2.0): '))
except ValueError:
    T = 2.0

af = 2*np.pi/T

tmax = 5*T                   # [s] duration time
dt = 0.1                   # [s] interval time
t_a = np.arange(0, tmax, dt)    # time after step stress
t_b = np.arange(-2.0,0,dt) # time before step stress
t = np.concatenate([t_b,t_a])   # whole time 
zeros = np.zeros(len(t_b))
e_a = np.array([eamp*np.sin(af*t) for t in t_a])
e = np.concatenate([zeros,e_a]) # whole strain

# solution of ODE（ここではそれをする必要はない）
s = E*e
el = e*l                        # [m] elongation

# scaling for figure
s = s/10**6                     # MPaスケール

fig = plt.figure(figsize=(8,5), tight_layout=True)
ax = fig.add_subplot(111)
ax.grid()
title_text = "spring (Hooke's elasticity): sinusoidal strain"
ax.set_title(title_text)
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
ax.set_xlim(-0.05,0.3)
ax.set_ylim(-5,5)

# for common
y_0 = [0, 0]
ax.plot([0, l/4],y_0, c='b')
ax.plot(0,0, 'ro', markersize='10')
ax.plot([0,l],[-2,-2], c='g')
ax.plot([0,0],[-1.8,-2.2], c='g')
ax.plot([l,l],[-1.8,-2.2], c='g')

var_text = r'$\epsilon_{{amp}}$ = {0:.2f}, $T$ = {1:.1f} s, $E$ = {2:.1f} MPa'.format(eamp,T,E/10**6)
ax.text(0.4, 0.9, var_text, transform=ax.transAxes)
eq_text = r'$\sigma$ = $E\epsilon$'
ax.text(0.4, 0.8, eq_text, transform=ax.transAxes)
ax.text(0.3, 0.25, '$l_0$', transform=ax.transAxes)
ax.text(0.6, 0.38, '$\sigma$ (output)', transform=ax.transAxes)

# for spring
rod, = ax.plot([],[], 'b', animated=True)
triangle, = ax.plot([],[], 'b', animated=True)
point, = ax.plot([],[], 'ro', markersize='10', animated=True)
# ここでは[],[]としているが、下で***.set_dataで実際の値を入れている

# for stress
stress, = ax.plot([],[], 'r', lw=2, animated=True)
arrow_p, = ax.plot([],[], 'r', marker=9, markersize='10', animated=True)
arrow_n, = ax.plot([],[], 'r', marker=8, markersize='10', animated=True)

time_template = '$t$ = %.1f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():               # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return rod, triangle, point, stress, arrow_p, arrow_n, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    x_rod = [3*l/4 + el[i], l + el[i]]
    rod.set_data(x_rod,y_0)
    x_tri = np.linspace(l/4, 3*l/4 + el[i],100)
    y_tri = w*((2/3)*np.arccos(np.cos(6*np.pi*(x_tri - l/4)/(el[i]+l/2)-np.pi/2+0.1))-1)
    triangle.set_data(x_tri,y_tri)
    point.set_data([l + el[i]],[0])
    a = 0.2    # 見かけ上の振幅
    x_stress = [l, l + a*s[i]/eamp]
    stress.set_data(x_stress,[-1,-1])
    if s[i] > 0:
        arrow_p.set_data([l + a*s[i]/eamp],[-1])
    else:
        arrow_n.set_data([l + a*s[i]/eamp],[-1])
    time_text.set_text(time_template % (t[i]))
    return rod, triangle, point, stress, arrow_p, arrow_n, time_text

f = np.arange(0, len(t))
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=f, 
                    init_func=init, blit=True, interval=frame_int, repeat=False)

savefile = "./gif/spring_sinuStrain_(epsilon={0:.2f},T={1:.1f}s,mod={2:.1f}MPa).gif".format(eamp,T,E/10**6)
ani.save(savefile, writer='pillow', fps=fps)

plt.show()