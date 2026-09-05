# response of 4-component model to time-dependent strain

import numpy as np
import matplotlib.pyplot as plt
import SLS2_time_dependent_strain as SLS2

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
        eta = float(input('viscosity [kPa s] (default = 100.0 kPa s): '))*10**3
    except ValueError:
        eta = 100*10**3               # [Pa s] 粘度
    try:
        Ep = float(input('probe modulus [MPa] (default = 5.0 MPa): '))*10**6
    except ValueError:
        Ep = 5*10**6                # [Pa] 弾性率

    # パラメータの計算
    insMod = E1 + E2             # [Pa] 瞬間弾性率
    infMod = E2                  # [Pa] 緩和弾性率
    relaxTime = eta/E1          # [s] 緩和時間
    k = insMod/infMod
    Er = (E2 * Ep) / (E2 + Ep)
    phi = (Ep/k + E2) / (E2 + Ep)
    tau_s = k * relaxTime * phi

    return E1, E2, eta, Ep, insMod, infMod, relaxTime, k, Er, phi, tau_s

def func_4comp_1(modulus, relaxTime, k, c, t1):
    # modulus must be Er
    # relaxTime must be tau_s
    # k must be 1/phi
    tim = np.linspace(0, t1, 400)
    strain = c*tim
    stress = (c*modulus)*(tim + relaxTime * (k - 1) * (1 - np.exp(-tim/relaxTime)))
    return tim, strain, stress

def func_4comp_2(modulus, relaxTime, c, t1, t2, s1):
    # modulus must be Er
    # relaxTime must be tau_s
    # k must be 1/phi
    tim = np.linspace(t1, t2, 400)
    strain = c*t1*np.ones(len(tim))
    stress = (c*modulus) * (t1 + (s1/(c*modulus) - t1) * np.exp(-(tim - t1)/relaxTime))
    return tim, strain, stress

def func_4comp_3(modulus, relaxTime, c, t1, t2, t3, s2, k):
    # modulus must be Er
    # relaxTime must be tau_s
    # k must be 1/phi
    tim = np.linspace(t2, t3, 400)
    strain = c*(t1 + t2 - tim)
    stress = (c*modulus)*(relaxTime * (1 - k) + t1 + t2 - tim + (s2/(c*modulus) - t1 - relaxTime * (1 - k)) * np.exp(-(tim-t2)/relaxTime))
    return tim, strain, stress

def func_4comp_4(relaxTime, t3, t4, s3):
    # modulus must be Er
    # relaxTime must be tau_s
    # k must be 1/phi
    tim = np.linspace(t3, t4, 400)
    strain = np.zeros(len(tim))
    stress = s3 * np.exp(-(tim-t3)/relaxTime)
    return tim, strain, stress

if __name__=='__main__':
# Parameters for 4-component model
    E1, E2, eta, Ep, insMod, infMod, relaxTime, k, Er, phi, tau_s = reqParams()
    param_text = r'($E_1$ = {0:.1f} MPa, $E_2$ = {1:.1f} MPa, $\eta$ = {2:.1f} kPa s)'.format(E1/10**6, E2/10**6, eta/10**3)
    probe_text = r'$E_p$ = {0:.1f} MPa'.format(Ep/10**6)
    res_text1 = r'$E_i$ = {0:.2f} MPa, $E_\infty$ = {1:.2f} MPa, $k$ = {2:.2f}, $\tau$ = {3:.2f} ms'.format(insMod/10**6, infMod/10**6, k, relaxTime*10**3)
    res_text2 = r'$E_r$ = {0:.2f} MPa, $\tau_s$ = {1:.2f} ms'.format(Er/10**6, tau_s*10**3)

# SLS2 model for comparison
    c, t1, tim1, strain1_sls, stress1_sls = SLS2.func_SLS2_1(infMod, relaxTime, k)
    s1_sls = stress1_sls[-1]    # stress at t1
    t2, tim2, strain2_sls, stress2_sls = SLS2.func_SLS2_2(infMod, relaxTime, c, t1, s1_sls)
    s2_sls = stress2_sls[-1]    # stress at t2
    t3, tim3, strain3_sls, stress3_sls = SLS2.func_SLS2_3(infMod, relaxTime, c, t1, t2, s2_sls, k)
    s3_sls = stress3_sls[-1]    # stress at t3
    t4, tim4, strain4_sls, stress4_sls = SLS2.func_SLS2_4(relaxTime, t3, s3_sls)

# 4-component model for time-dependent strain
    tim1, strain1, stress1 = func_4comp_1(Er, tau_s, 1/phi, c, t1)
    s1 = stress1[-1]    # stress at t1
    tim2, strain2, stress2 = func_4comp_2(Er, tau_s, c, t1, t2, s1)
    s2 = stress2[-1]    # stress at t2
    tim3, strain3, stress3 = func_4comp_3(Er, tau_s, c, t1, t2, t3, s2, 1/phi)
    s3 = stress3[-1]    # stress at t3
    tim4, strain4, stress4 = func_4comp_4(tau_s, t3, t4, s3)

    tim0 = np.linspace(-t1, 0, 400)
    stress0 = np.zeros(len(tim0))
    strain0 = np.zeros(len(tim0))

    tim = np.concatenate([tim0,tim1,tim2,tim3,tim4])/10**(-3)   # rescale to ms
    strain = np.concatenate([strain0,strain1,strain2,strain3,strain4])
    stress = np.concatenate([stress0,stress1,stress2,stress3,stress4])/10**6   # rescale to MPa
    stress_sls = np.concatenate([stress0,stress1_sls,stress2_sls,stress3_sls,stress4_sls])/10**6   # rescale to MPa

    strain_probe = stress / (Ep/10**6)   # rescale to MPa
    strain_mat = strain - strain_probe

    try:
        select = int(input('Selection (strain&stress vs time: 0, stress vs strain: 1): '))
    except ValueError:
        select = 0

    if select == 0:
        fig = plt.figure(figsize=(8,10), tight_layout=True)
        ax1 = fig.add_subplot(211)
        ax1.set_title('4-component model for time-dependent strain '+param_text)
        ax1.set_xlabel(r'$t$ /ms')
        ax1.set_ylabel(r'$\epsilon$ /')
        ax1.set_xlim(-t1*10**3, t4*10**3)
        ymin = np.maximum(np.abs(1.4*np.min(strain_probe)), 0.1*np.max(strain))
        ax1.set_ylim(-ymin, 1.4*np.max(strain))
        ax1.grid()
        ax1.set_axisbelow(True)
        ax1.text(0.05, 0.95, probe_text, transform=ax1.transAxes, verticalalignment='top')
        ax1.text(0.05, 0.90, res_text1, transform=ax1.transAxes, verticalalignment='top')
        ax1.text(0.05, 0.85, res_text2, transform=ax1.transAxes, verticalalignment='top')
        ax1.plot(tim, strain, c='b', lw=1, ls='--', label='strain (input)')
        ax1.plot(tim, strain_probe, c='g', lw=2, ls='--', label='strain (probe)')
        ax1.plot(tim, strain_mat, c='cyan', lw=2, ls='-', label='strain (material)')
        ax1.legend(loc='upper right')

        ax2 = fig.add_subplot(212)
        ax2.set_xlabel(r'$t$ /ms')
        ax2.set_ylabel(r'$\sigma$ /MPa')
        ax2.set_xlim(-t1*10**3, t4*10**3)
        ax2.set_ylim(1.2*np.min(stress_sls), 1.2*np.max(stress_sls))
        ax2.grid()
        ax2.set_axisbelow(True)
        ax2.plot(tim, stress, c='r', lw=2, label='Response to time-dependent strain')
        ax2.plot(tim, stress_sls, c='black', lw=1, ls='--', label='Response to time-dependent strain (SLS2)')
        ax2.legend(loc='upper right')
        savefile = './png/4-component_time-dependent_strain_(Ep={0:.1f}MPa).png'.format(Ep/10**6)

    elif select == 1:
        fig = plt.figure(figsize=(8,5), tight_layout=True)
        ax1 = fig.add_subplot(111)
        ax1.set_title('4-component model for time-dependent strain '+param_text)
        ax1.set_xlabel(r'$\epsilon$ /')
        ax1.set_ylabel(r'$\sigma$ /MPa')
        ax1.set_xlim(-0.05*np.max(strain), 1.2*np.max(strain))
        ax1.set_ylim(1.2*np.min(stress_sls), 1.2*np.max(stress_sls))
#        ax1.set_ylim(1.2*np.min(stress), 1.2*np.max(stress))
        ax1.grid()
        ax1.set_axisbelow(True)
        ax1.text(0.05, 0.95, probe_text, transform=ax1.transAxes, verticalalignment='top')
        ax1.text(0.05, 0.90, res_text1, transform=ax1.transAxes, verticalalignment='top')
        ax1.text(0.05, 0.85, res_text2, transform=ax1.transAxes, verticalalignment='top')
        ax1.plot(strain, stress, c='r', lw=1, ls='--', label='stress-strain curve')
        ax1.plot(strain_mat, stress, c='violet', lw=2, label='stress-strain curve (material)')
        ax1.plot(strain_probe, stress, c='b', lw=1, ls='--', label='stress-strain curve (probe)')
        ax1.plot(strain, stress_sls, c='black', lw=1, ls='--', label='stress-strain curve (SLS2)')
        ax1.legend(loc='upper right')
        savefile = './png/4-component_stress-strain_curve_(Ep={0:.1f}MPa).png'.format(Ep/10**6)

    fig.savefig(savefile, dpi=300)
    plt.show()