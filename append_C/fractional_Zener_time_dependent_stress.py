# response of SLS1 model to time-dependent stress

import numpy as np
import matplotlib.pyplot as plt
from differintP.functions import MittagLeffler

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
        E = float(input('Enter modulus value of spring-pot (MPa) (default = 0.8 MPa): '))*10**6
    except ValueError:
        E = 8*10**5
    try:
        eta = float(input('Enter viscosity value of spring-pot (kPa s) (default = 100 kPa s): '))*10**3
    except ValueError:
        eta = 10**5

    return E1, E2, E, eta

def timeAxes():
    try:
        c = float(input('Enter c of stress = c*time (default = 0.01 MPa/ms): '))
    except ValueError:
        c = 10**7
    try:
        t1 = float(input('Enter t1 (0<=t<t1) (ms) (default = 40 ms): '))*10**(-3)
    except ValueError:
        t1 = 40*10**(-3)
    try:
        dt1 = float(input('Enter dt = t2 - t1 (t1<=t<t2) (default = 40 ms): '))*10**(-3)
    except ValueError:
        dt1 = 40*10**(-3)
    try:
        dt2 = float(input('Enter dt = t3 - t2 (t2<=t<t3) (default = 40 ms): '))*10**(-3)  
    except ValueError:
        dt2 = 40*10**(-3)
    try:
        dt3 = float(input('Enter dt = t4 - t3 (t3<=t<t4) (default = 80 ms): '))*10**(-3)  
    except ValueError:
        dt3 = 80*10**(-3)
    t2 = t1 + dt1
    t3 = t2 + dt2
    t4 = t3 + dt3   
    tim1 = np.linspace(0, t1, 400)
    tim2 = np.linspace(t1, t2, 400)
    tim3 = np.linspace(t2, t3, 400)
    tim4 = np.linspace(t3, t4, 400) 
    return c, t1, t2, t3, t4, tim1, tim2, tim3, tim4

def func_fracZ_1(c, tim, modulus, tau, k, nu):
    # modulus must be infMod
    # tau must be tau_prime
    # tim must be tim1
    x0 = -(1/k)*(tim/tau)**nu
    stress = c*tim
    strain = (c/modulus)*tim*(1 - (1 - 1/k)*MittagLeffler(nu, 2, x0, num_terms=300))
    strain_sls2 = (c/modulus)*(tim - k*tau*(1 - 1/k)*(1 - np.exp(-tim/(k*tau))))
    return stress, strain, strain_sls2

def func_fracZ_2(c, tim, modulus, tau, t1, k, nu):
    # modulus must be infMod
    # tau must be tau_prime
    # tim must be tim2
    x0 = -(1/k)*(tim/tau)**nu
    x1 = -(1/k)*((tim - t1)/tau)**nu
    stress = c*t1*np.ones(len(tim))
    strain = (c/modulus) * (t1 - (1 - 1/k)*(tim*MittagLeffler(nu, 2, x0, num_terms=300) - (tim - t1)*MittagLeffler(nu, 2, x1, num_terms=300)))
    strain_sls2 = (c/modulus)*(t1 - k*tau*(1 - 1/k)*(1 - np.exp(-t1/(k*tau)))*np.exp(-(tim - t1)/(k*tau)))
    return stress, strain, strain_sls2

def func_fracZ_3(c, tim, modulus, tau, t1, t2, k, nu):
    # modulus must be infMod
    # tau must be tau_prime
    # tim must be tim3
    x0 = -(1/k)*(tim/tau)**nu
    x1 = -(1/k)*((tim - t1)/tau)**nu
    x2 = -(1/k)*((tim - t2)/tau)**nu
    stress = c*(t1 + t2 - tim)
    strain = (c/modulus) * (t1 + t2 - tim - (1 - 1/k)*(tim*MittagLeffler(nu, 2, x0, num_terms=300) - (tim - t1)*MittagLeffler(nu, 2, x1, num_terms=300) - (tim - t2)*MittagLeffler(nu, 2, x2, num_terms=300)))
    strain_sls2 = (c/modulus)*((k*tau*(1 - 1/k) + t1 + t2 - tim) - k*tau*(1 - 1/k)*(1 + (1 - np.exp(-t1/(k*tau)))*np.exp(-(t2 - t1)/(k*tau)))*np.exp(-(tim - t2)/(k*tau)))
    return stress, strain, strain_sls2

def func_fracZ_4(c, tim, modulus, tau, t1, t2, t3, k, nu):
    # modulus must be infMod
    # tau must be tau_prime
    # tim must be tim4
    x0 = -(1/k)*(tim/tau)**nu
    x1 = -(1/k)*((tim - t1)/tau)**nu
    x2 = -(1/k)*((tim - t2)/tau)**nu
    x3 = -(1/k)*((tim - t3)/tau)**nu
    stress = np.zeros(len(tim))
    strain = -(c/modulus) * (1 - 1/k)*(tim*MittagLeffler(nu, 2, x0, num_terms=300) - (tim - t1)*MittagLeffler(nu, 2, x1, num_terms=300) - (tim - t2)*MittagLeffler(nu, 2, x2, num_terms=300) + (tim - t3)*MittagLeffler(nu, 2, x3, num_terms=300))
    strain_sls2 = (c/modulus)*k*tau*(1 - 1/k)*(1 - (1 + (1 - np.exp(-t1/(k*tau)))*np.exp(-(t2 - t1)/(k*tau)))*np.exp(-(t3 - t2)/(k*tau)))*np.exp(-(tim - t3)/(k*tau))
    return stress, strain, strain_sls2

cmap = plt.get_cmap('winter')
cmap2 = plt.get_cmap('autumn')

if __name__=='__main__':
    E1, E2, E, eta = reqParams()
    insMod = E1 + E2             # [Pa] 瞬間弾性率
    infMod = E2                  # [Pa] 緩和弾性率
    k = insMod/infMod
    tau = eta/E                  # [s] 緩和時間    
    kappa = E/insMod

    param_text = """
    ($E_1$ = {0:.1f} MPa, $E_2$ = {1:.1f} MPa, $E$ = {2:.1f} MPa, $\eta$ = {3:.1f} kPa s)""".format(E1/10**6, E2/10**6, E/10**6, eta/10**3)
    res_text = r'$E_i$ = {0:.2f} MPa, $E_\infty$ = {1:.2f} MPa, $k$ = {2:.2f}, $\tau$ = {3:.2f} ms'.format(insMod/10**6, infMod/10**6, k, tau*10**3)

    c, t1, t2, t3, t4, tim1, tim2, tim3, tim4 = timeAxes()

    tim0 = np.linspace(-t1, 0, 400)
    stress0 = np.zeros(len(tim0))
    strain0 = np.zeros(len(tim0))

    tim = np.concatenate([tim0,tim1,tim2,tim3,tim4])/10**(-3)   # rescale to ms
    nu_arr = np.array([0.2, 0.4, 0.6, 0.8, 1.0])

    stress_array = np.zeros((len(nu_arr), len(tim)))  # 階数ごとの応力を格納する配列
    strain_array = np.zeros((len(nu_arr), len(tim)))  # 階数ごとの歪みを格納する配列
    strain_sls2_array = np.zeros((len(nu_arr), len(tim)))  # 階数ごとの歪みを格納する配列
    tau_prime_arr = np.zeros(len(nu_arr))  # 階数ごとの緩和時間を格納する配列

    for i in range(len(nu_arr)):
        nu = nu_arr[i]
        tau_prime = kappa**(1/nu)*tau
        tau_prime_arr[i] = tau_prime
        stress1, strain1, strain_sls2_1 = func_fracZ_1(c, tim1, infMod, tau_prime, k, nu)
        stress2, strain2, strain_sls2_2 = func_fracZ_2(c, tim2, infMod, tau_prime, t1, k, nu)
        stress3, strain3, strain_sls2_3 = func_fracZ_3(c, tim3, infMod, tau_prime, t1, t2, k, nu)
        stress4, strain4, strain_sls2_4 = func_fracZ_4(c, tim4, infMod, tau_prime, t1, t2, t3, k, nu)
        stress = np.concatenate([stress0,stress1,stress2,stress3,stress4])/10**6   # rescale to MPa
        strain = np.concatenate([strain0,strain1,strain2,strain3,strain4])
        strain_sls2 = np.concatenate([strain0,strain_sls2_1,strain_sls2_2,strain_sls2_3,strain_sls2_4])
        stress_array[i] = stress
        strain_array[i] = strain
        strain_sls2_array[i] = strain_sls2

    try:
        select = int(input('Selection (strain&stress vs time: 0, stress vs strain: 1): '))
    except ValueError:
        select = 0

    if select == 0:
        title = 'fractinal Zener model for time-dependent stress '+param_text
        x_label = r'$t$ /ms'
        y1_label = r'$\sigma$ /MPa'
        y2_label = r'$\epsilon$ /'
        x_lim = [np.min(tim),np.max(tim)]
        y1_lim = [-0.1*np.max(stress_array), 1.2*np.max(stress_array)]
        y2_lim = [-0.1*np.nanmax(strain_array), 1.2*np.nanmax(strain_array)]
        legend_loc='upper right'
        savefile = './png/fracZ_time-dependent_stress_(tau={0:.1f}ms).png'.format(tau*10**3)

    if select == 1:
        title = 'Stress-strain curve for fractional Zener model '+param_text
        x_label = r'$\epsilon$ /'
        y_label = r'$\sigma$ /MPa'
        x_lim = [-0.1*np.nanmax(strain_array), 1.5*np.nanmax(strain_array)]
        y_lim = [-0.1*np.max(stress_array), 1.2*np.max(stress_array)]
        legend_loc='upper right'
        savefile = './png/fracZ_stress-strain_curve_(tau={0:.1f}ms).png'.format(tau*10**3)

    if select == 0:
        fig = plt.figure(figsize=(8,10), tight_layout=True)
        ax1 = fig.add_subplot(211)
        ax1.set_title(title)
        ax1.set_xlabel(x_label)
        ax1.set_ylabel(y1_label)
        ax1.set_xlim(x_lim[0], x_lim[1])
        ax1.set_ylim(y1_lim[0], y1_lim[1])
        ax1.grid()
        ax1.set_axisbelow(True)
        ax1.text(0.05, 0.95, res_text, transform=ax1.transAxes, verticalalignment='top')
        ax1.plot(tim, stress_array[0], c='r', lw=2, label='Time-dependent stress')
        ax1.legend(loc='upper right')

        ax2 = fig.add_subplot(212)
        ax2.set_xlabel(x_label)
        ax2.set_ylabel(y2_label)
        ax2.set_xlim(x_lim[0], x_lim[1])
        ax2.set_ylim(y2_lim[0], y2_lim[1])
        ax2.grid()
        ax2.set_axisbelow(True)
        for i in range(len(nu_arr)):
            ax2.plot(tim, strain_array[i], color=cmap(i/len(nu_arr)), label=r'$\nu$ = {0:.1f}, $\tau^\prime$ = {1:.1f} $\mu$s'.format(nu_arr[i], tau_prime_arr[i]*10**6))
        # 関数にtau_primeを渡しているので、tauとtau_primeが一致するnu=1.0のときのSLS IIモデルの応答をプロットする
        ax2.plot(tim, strain_sls2_array[4], c='r', lw=1, ls='--', label='SLS II')
        ax2.legend(loc='upper right')

    elif select == 1:
        fig = plt.figure(figsize=(8,5), tight_layout=True)
        ax1 = fig.add_subplot(111)
        ax1.set_title(title)
        ax1.set_xlabel(x_label)
        ax1.set_ylabel(y_label)
        ax1.set_xlim(x_lim[0], x_lim[1])
        ax1.set_ylim(y_lim[0], y_lim[1])
        ax1.grid()
        ax1.set_axisbelow(True)
        ax1.text(0.05, 0.95, res_text, transform=ax1.transAxes, verticalalignment='top')
        for i in range(len(nu_arr)):
            ax1.plot(strain_array[i], stress_array[i], color=cmap2(i/len(nu_arr)), label=r'$\nu$ = {0:.1f}, $\tau^\prime$ = {1:.1f} $\mu$s'.format(nu_arr[i], tau_prime_arr[i]*10**6))
        ax1.plot(strain_sls2_array[4], stress_array[0], c='green', lw=1, ls='--', label='SLS II')
        ax1.legend(loc='upper right')
        savefile = './png/fracZ_stress-strain_curve_(tau={0:.1f}ms).png'.format(tau*10**3)

    fig.savefig(savefile, dpi=300)
    plt.show()