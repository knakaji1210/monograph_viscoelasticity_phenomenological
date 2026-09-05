# response of SLS1 model to time-dependent stress

import numpy as np
import matplotlib.pyplot as plt

def reqParams():
    # 変数の設定
    try:
        E1 = float(input('modulus 1 [MPa] (default = 10.0 MPa): '))*10**6
    except ValueError:
        E1 = 10**7                  # [Pa] 弾性率
    try:
        E2 = float(input('modulus 2 [MPa] (default = 1.0 MPa): '))*10**6
    except ValueError:
        E2 = 10**6                # [Pa] 弾性率
    try:
        eta = float(input('viscosity [kPa s] (default = 10.0 kPa s): '))*10**3
    except ValueError:
        eta = 10*10**3               # [Pa s] 粘度

    # パラメータの計算
    insMod = E1                 # [Pa] 瞬間弾性率
    infMod = E1*E2/(E1+E2)      # [Pa] 緩和弾性率
    retardTime = eta/E2          # [s] 緩和時間
    k = insMod/infMod

    return E1, E2, eta, insMod, infMod, retardTime, k

def func_SLS1_1(modulus, retardTime, k):
    # modulus must be infMod
    try:
        c = float(input('Enter c of stress = c*time (default = 0.01 MPa/ms): '))
    except ValueError:
        c = 10**7
    try:
        t1 = float(input('Enter t1 (0<=t<t1) (ms) (default = 40 ms): '))*10**(-3)
    except ValueError:
        t1 = 40*10**(-3)
    tim = np.linspace(0, t1, 400)
    stress = c*tim
    strain = (c/modulus)*(tim - retardTime * (1 - 1/k) * (1 - np.exp(-tim/retardTime)))
    return c, t1, tim, stress, strain

def func_SLS1_2(modulus, retardTime, c, t1, e1):
    # modulus must be infMod
    try:
        dt = float(input('Enter dt = t2 - t1 (t1<=t<t2) (default = 60 ms): '))*10**(-3)
    except ValueError:
        dt = 60*10**(-3)
    t2 = t1 + dt
    tim = np.linspace(t1, t2, 400)
    stress = c*t1*np.ones(len(tim))
    strain = (c/modulus) * (t1 + (e1/(c/modulus) - t1) * np.exp(-(tim - t1)/retardTime))
    return t2, tim, stress, strain

def func_SLS1_3(modulus, retardTime, c, t1, t2, e2, k):
    # modulus must be infMod
    try:
        dt = float(input('Enter dt = t3 - t2 (t2<=t<t3) (default = 40 ms): '))*10**(-3)  
    except ValueError:
        dt = 40*10**(-3)
    t3 = t2 + dt
    tim = np.linspace(t2, t3, 400)
    stress = c*(t1 + t2 - tim)
    strain = (c/modulus)*(retardTime * (1 - 1/k) + t1 + t2 - tim + (e2/(c/modulus) - t1 - retardTime * (1 - 1/k)) * np.exp(-(tim-t2)/retardTime))
    return t3, tim, stress, strain

def func_SLS1_4(retardTime, t3, e3):
    # modulus must be infMod
    try:
        dt = float(input('Enter dt = t4 - t3 (t3<=t<t4) (default = 40 ms): '))*10**(-3)  
    except ValueError:
        dt = 40*10**(-3)
    t4 = t3 + dt
    tim = np.linspace(t3, t4, 400)
    stress = np.zeros(len(tim))
    strain = e3 * np.exp(-(tim-t3)/retardTime)
    return t4, tim, stress, strain

if __name__=='__main__':
    E1, E2, eta, insMod, infMod, retardTime, k = reqParams()
    param_text = r'($E_1$ = {0:.1f} MPa, $E_2$ = {1:.1f} MPa, $\eta$ = {2:.1f} kPa s)'.format(E1/10**6, E2/10**6, eta/10**3)
    res_text = r'$E_i$ = {0:.2f} MPa, $E_\infty$ = {1:.2f} MPa, $k$ = {2:.2f}, $\tau$ = {3:.2f} ms'.format(insMod/10**6, infMod/10**6, k, retardTime*10**3)
    c, t1, tim1, stress1, strain1 = func_SLS1_1(infMod, retardTime, k)
    e1 = strain1[-1]    # strain at t1
    t2, tim2, stress2, strain2 = func_SLS1_2(infMod, retardTime, c, t1, e1)
    e2 = strain2[-1]    # strain at t2
    t3, tim3, stress3, strain3 = func_SLS1_3(infMod, retardTime, c, t1, t2, e2, k)
    e3 = strain3[-1]    # strain at t3
    t4, tim4, stress4, strain4 = func_SLS1_4(retardTime, t3, e3)

    tim0 = np.linspace(-t1, 0, 400)
    stress0 = np.zeros(len(tim0))
    strain0 = np.zeros(len(tim0))

    tim = np.concatenate([tim0,tim1,tim2,tim3,tim4])/10**(-3)   # rescale to ms
    stress = np.concatenate([stress0,stress1,stress2,stress3,stress4])/10**6   # rescale to MPa
    strain = np.concatenate([strain0,strain1,strain2,strain3,strain4])

    try:
        select = int(input('Selection (strain&stress vs time: 0, stress vs strain: 1): '))
    except ValueError:
        select = 0

    if select == 0:
        fig = plt.figure(figsize=(8,10), tight_layout=True)
        ax1 = fig.add_subplot(211)
        ax1.set_title('SLS1 model for time-dependent strain '+param_text)
        ax1.set_xlabel(r'$t$ /ms')
        ax1.set_ylabel(r'$\sigma$ /MPa')
        ax1.set_xlim(-t1*10**3, t4*10**3)
        ax1.set_ylim(-0.05*np.max(stress), 1.2*np.max(stress))
        ax1.grid()
        ax1.set_axisbelow(True)
        ax1.text(0.05, 0.95, res_text, transform=ax1.transAxes, verticalalignment='top')
        ax1.plot(tim, stress, c='r', lw=2, label='Time-dependent stress')
        ax1.legend(loc='upper right')

        ax2 = fig.add_subplot(212)
        ax2.set_xlabel(r'$t$ /ms')
        ax2.set_ylabel(r'$\epsilon$ /')
        ax2.set_xlim(-t1*10**3, t4*10**3)
        ax2.set_ylim(-0.05*np.max(strain), 1.2*np.max(strain))
        ax2.grid()
        ax2.set_axisbelow(True)
        ax2.plot(tim, strain, c='b', lw=2, label='Response to time-dependent stress')
        ax2.legend(loc='upper right')
        savefile = './png/SLS1_time-dependent_strain_(tau={0:.1f}ms).png'.format(retardTime*10**3)

    elif select == 1:
        fig = plt.figure(figsize=(8,5), tight_layout=True)
        ax1 = fig.add_subplot(111)
        ax1.set_title('SLS1 model for time-dependent strain '+param_text)
        ax1.set_xlabel(r'$\epsilon$ /')
        ax1.set_ylabel(r'$\sigma$ /MPa')
        ax1.set_xlim(-0.05*np.max(strain), 1.2*np.max(strain))
        ax1.set_ylim(-0.05*np.max(stress), 1.2*np.max(stress))
        ax1.grid()
        ax1.set_axisbelow(True)
        ax1.text(0.05, 0.95, res_text, transform=ax1.transAxes, verticalalignment='top')
        ax1.plot(strain, stress, c='r', lw=2, label='stress-strain curve')
        ax1.legend(loc='upper right')
        savefile = './png/SLS1_stress-strain_curve_(tau={0:.1f}ms).png'.format(retardTime*10**3)

    fig.savefig(savefile, dpi=300)
    plt.show()