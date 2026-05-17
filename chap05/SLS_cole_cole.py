# Cole-Cole plot (SLS I & SLS II)
import numpy as np
import matplotlib.pyplot as plt

def reqParams(model):
    xi_min, xi_max = -2, 2
    xi = np.logspace(xi_min, xi_max, 200)
    # 変数の設定
    try:
        E1 = float(input('modulus 1 [MPa] (default = 1.0 MPa): '))*10**6
    except ValueError:
        E1 = 10**6                  # [Pa] 弾性率
    try:
        E2 = float(input('modulus 2 [MPa] (default = 0.2 MPa): '))*10**6
    except ValueError:
        E2 = 2*10**5                # [Pa] 弾性率
    if model == 0:
        # 複素弾性率の計算
        insMod = E1                 # [Pa] 瞬間弾性率
        infMod = E1*E2/(E1+E2)      # [Pa] 緩和弾性率
        k = insMod/infMod    
        numer = insMod*(1 + xi*(2j/2))
        denom = k + xi*(2j/2)
        comMod = numer/denom
        x_label = r'$J^{\prime}$ / $J_i$'
        y_label = r'$J^{{\prime\prime}}$ / $J_i$'
        modeltext = r'SLS I model '
        savefile = './png/SLS1_cole_cole.png'

    if model == 1:
        # 複素弾性率の計算
        insMod = E1+E2              # [Pa] 瞬間弾性率
        infMod = E2                 # [Pa] 緩和弾性率
        k = insMod/infMod              # [s] 緩和時間    
        numer = insMod*(1/k + xi*(2j/2))
        denom = 1 + xi*(2j/2)
        comMod = numer/denom
        x_label = r'$E^{\prime}$ / $E_i$'
        y_label = r'$E^{{\prime\prime}}$ / $E_i$'
        modeltext = r'SLS II model '
        savefile = './png/SLS2_cole_cole.png'

    return k, insMod, infMod, comMod, x_label, y_label, modeltext, savefile

if __name__=='__main__':
    try:
        model = int(input('Selection (SLS I : 0, SLS II: 1): '))
    except ValueError:
        model = 0    
    k, insMod, infMod, comMod, x_label, y_label, modeltext, savefile = reqParams(model)    
    param_text = '($k$ = {0})'.format(k)

    if model == 0:
        comComp = 1 / comMod
        strComp = comComp.real
        losComp = -comComp.imag
        x = strComp * insMod
        y = losComp * insMod

    elif model == 1:
        strMod = comMod.real
        losMod = comMod.imag
        x = strMod / insMod
        y = losMod / insMod

    fig = plt.figure(figsize=(8,8), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.set_title(modeltext+param_text)
    ax.set_xlim(-0.05*np.max(x), 1.05*np.max(x))
    ax.set_ylim(-0.05*np.max(x), 1.05*np.max(x))
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.scatter(x, y)
    ax.grid()
    ax.set_axisbelow(True)

fig.savefig(savefile, dpi=300)

plt.show()