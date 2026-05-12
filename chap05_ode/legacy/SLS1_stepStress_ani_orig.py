# ordinary differential equation of Voigt model (step stress) with animation

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.animation import FuncAnimation

'''
テキストの式(4.11)をベースに組み立てる
'''

def SLS1_stepStress(e, t, s, infMod, tau):
# e: strain, s: stress, infMod: relaxation modulus, tau: retardation time
# E1: insMod, E1*E2/(E1+E2): infMod
# ここでは下でargsとしてs=s0を入れてステップ応力を実現
    dedt = (s/infMod - e)/tau   # (4.11)
    return dedt

# variables
try:
    E1 = float(input('modulus 1 [MPa] (default = 1.0 MPa): '))*10**6
except ValueError:
    E1 = 10**6                  # [Pa] modulus
try:
    E2 = float(input('modulus 2 [MPa] (default = 0.5 MPa): '))*10**6
except ValueError:
    E2 = 5*10**5                # [Pa] modulus
try:
    eta = float(input('viscosity [kPa s] (default = 1000.0 kPa s): '))*10**3
except ValueError:
    eta = 10**6             # [Pa s] viscosity
tau = eta/E2                  # [s] retardation time
insMod = E1                   # [Pa] instantaneous modulus
infMod = E1*E2/(E1+E2)        # [Pa] equilibrium modulus
k = insMod/infMod
l = 0.1                     # [m] equilibrium length
w = 0.5                     # ratio of dashpot width

# initial condition
try:
    s0 = float(input('step stress [MPa] (default = 0.1 MPa): '))*10**6
except ValueError:
    s0 = 10**5              # [Pa] step stress
e0 = s0/E1                  # [] initial strain

tmax = 10                   # [s] duration time
dt = 0.05                   # [s] interval time
t_a = np.arange(0, tmax, dt)    # time after step stress
t_b = np.arange(-2.0,0,dt) # time before step stress
t = np.concatenate([t_b,t_a])   # whole time 
zeros = np.zeros(len(t_b))
ones = np.ones(len(t_a))
s = np.concatenate([zeros,ones*s0])

# solution of ODE
sol = odeint(SLS1_stepStress, e0, t_a, args=(s0,infMod,tau))
e = np.concatenate([zeros,sol[:,0]])        # [] strain
e1 = s/E1                                   # [] strain of spring1
e2 = e - e1                                 # [] strain of Voigt component
s1 = E2*e2                                  # [MPa] stress on sping in Voigt component
de2dt = np.array([0.0]+[(e2[k+1]-e2[k])/(t[k+1]-t[k]) for k in range(len(e)-1)])     # 簡易的なde/dt (Voigt部分)
s2 = eta*de2dt                              # [MPa] stress on dashpot
sv = s1 + s2                                # [MPa] Voigt部分の応力確認用

el_a = e*2*l      # 直列なのでl0=2lとなっている
el = e2*l         # Voigt要素部分

# scaling for figure
s = s/10**6                   # MPaスケール
s1 = s1/10**6                 # MPaスケール
s2 = s2/10**6                 # MPaスケール
sv = sv/10**6                 # MPaスケール

fig = plt.figure(figsize=(8,5), tight_layout=True)
ax = fig.add_subplot(111, xlabel='$t$ /s')
ax.grid()
title_text = "SLS1 model: step stress"
ax.set_title(title_text)
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
ax.set_xlim(-0.05,0.3)
ax.set_ylim(-5,5)
# for common
y_0 = [0, 0]
y_1 = [1, 1]
y_2 = [-1, -1]
ax.plot([0, l/4],y_0, c='b')
ax.plot([l/4, 9*l/4],[-2,-2], c='g')
ax.plot([l/4, l/4],[-1.8,-2.2], c='g')
ax.plot([9*l/4, 9*l/4],[-1.8,-2.2], c='g')
ax.plot([l/4,l/4],[1,-1], c='b')
ax.plot([l/4, (0.08+1/4)*l],y_1, c='b')
ax.plot([l/4, l/2],y_2, c='b')
ax.plot(0,0,'ro', markersize='10')
# for dashpot
x_d1 = [(0.08+1/4)*l, (0.92+1/4)*l]
y_d1 = [w+1, w+1]
y_d2 = [-w+1, -w+1]
ax.plot(x_d1,y_d1, c='b')
ax.plot(x_d1,y_d2, c='b')
ax.plot([(0.08+1/4)*l,(0.08+1/4)*l],[w+1,-w+1], c='b')
rect = patches.Rectangle(xy=((0.08+1/4)*l, -w+1), width=0.83*l, height=2*w, facecolor='y')
ax.add_patch(rect)

var_text = r'$\sigma_0$ = {0:.2f} MPa, $E_1$ = {1:.1f} MPa, $E_2$ = {2:.1f} MPa, $\eta$ = {3:.1f} kPa s'.format(s0/10**6,E1/10**6,E2/10**6,eta/10**3)
ax.text(0.4, 0.9, var_text, transform=ax.transAxes)
eq_text = r'd$\epsilon$/d$t$ = ($\sigma_0$/$E_{{{\infty}}}$ - $\epsilon$)/$\tau$'
ax.text(0.4, 0.8, eq_text, transform=ax.transAxes)
tau = eta/E2
res_text = r'$\tau$ = {0:.1f} s'.format(tau)
ax.text(0.4, 0.7, res_text, transform=ax.transAxes)
ax.text(0.5, 0.25, '$l_0$', transform=ax.transAxes)

# for common
bar, = ax.plot([],[], 'b', animated=True)
rod, = ax.plot([],[], 'b', animated=True)
point, = ax.plot([],[], 'ro', markersize='10', zorder=10, animated=True)
# for spring2
rod_sp2, = ax.plot([],[], 'b', animated=True)
triangle2, = ax.plot([],[], 'b', animated=True)
# for dashpot
rod_da, = ax.plot([],[], 'b', animated=True)
damper, = ax.plot([],[], 'b', lw=4, animated=True)
# for spring1
rod_sp1, = ax.plot([],[], 'b', animated=True)
triangle1, = ax.plot([],[], 'b', animated=True)

# ここでは[],[]としているが、下で***.set_dataで実際の値を入れている

time_template = '$t$ = %.1f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():               # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return bar, rod, point, rod_sp1, rod_sp2, triangle1, triangle2, rod_da, damper, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    # for common
    x_o = l + el[i]
    x_c = l/4 + x_o
    x_p = 9*l/4+el_a[i]
    bar.set_data([x_c,x_c],[1,-1])
    x_rod = [x_c, l/4+x_c]
    rod.set_data(x_rod,y_0)
    # for spring2
    x_rod_sp2 = [x_o, x_c]
    rod_sp2.set_data(x_rod_sp2,y_2)
    x_tri2 = np.linspace(l/2, x_o,100)
    y_tri2 = w*((2/3)*np.arccos(np.cos(6*np.pi*(x_tri2 - l/2)/(el[i]+l/2)-np.pi/2+0.1))-3)
    triangle2.set_data(x_tri2,y_tri2)
    # for dashpot
    x_rod_da = [x_o-l/4, x_c]
    x_damp = x_o-l/4
    y_damp = 0.7*w
    x_damper = [x_damp, x_damp]
    y_damper = [y_damp+1, -y_damp+1]
    rod_da.set_data(x_rod_da,y_1)
    damper.set_data(x_damper,y_damper)
    # for spring1
    x_rod_sp1 = [x_p-l/4 , x_p]
    rod_sp1.set_data(x_rod_sp1,y_0)
    x_tri1 = np.linspace(x_c+l/4, x_p-l/4,100)
    y_tri1 = w*((2/3)*np.arccos(np.cos(6*np.pi*(x_tri1 - x_c - l/4)/(el_a[i]-el[i]+l/2)-np.pi/2+0.1))-1)
    triangle1.set_data(x_tri1,y_tri1)
    point.set_data([x_p],[0])    
    time_text.set_text(time_template % (t[i]))
    return bar, rod, point, rod_sp1, rod_sp2, triangle1, triangle2, rod_da, damper, time_text

'''
y_triの中の重要部分は
x_tri1 = np.linspace(a, b,100)
のとき
(xtri - a)/(b - a)
になる 
'''

f = np.arange(0, len(t))
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=f, 
                    init_func=init, blit=True, interval=frame_int, repeat=False)

savefile = "./gif/SLS1_stepStress_(sigma={0:.2f}MPa,E1={1:.1f}MPa,E2={2:.1f}MPa,eta={3:.1f}kPas).gif".format(s0/10**6,E1/10**6,E2/10**6,eta/10**3)
ani.save(savefile, writer='pillow', fps=fps)

plt.show()