# Plemelj formula
import numpy as np
import matplotlib.pyplot as plt

def calc_Plemelj_real(a, x):
    y = a/(a**2 + x**2)
    return y

def calc_Plemelj_imaginary(a, x):
    y = x/(a**2 + x**2)
    return y

if __name__=='__main__':
    select_text = 'Selection (Real part: 0, Imaginary part: 1): '  
    try:
        select = int(input(select_text))
    except ValueError:
        select = 0

    a = np.logspace(-2, 0, 3)
    x = np.linspace(-3, 3, 1000)

    fig = plt.figure(tight_layout=True)
    ax = fig.add_subplot(111)
    ax.grid()
    ax.set_xlabel(r'$\omega$')

    if select == 0:
        for i in range(len(a)):
            title_text = 'Real part of Plemelj formula'
            savefile = './png/Plemelj_real.png'
            y = calc_Plemelj_real(a[i], x)
            ax.set_ylim(0,20)
            ax.plot(x, y, label=r'$\alpha$ = {0:.2f}'.format(a[i]))

    elif select == 1:
        for i in range(len(a)):
            title_text = 'Imaginary part of Plemelj formula'
            savefile = './png/Plemelj_imaginary.png'
            y = calc_Plemelj_imaginary(a[i], x)
            ax.set_ylim(-10,10)
            ax.plot(x, y, label=r'$\alpha$ = {0:.2f}'.format(a[i]))

    ax.set_title(title_text)
    ax.legend(loc='upper left')

    fig.savefig('./png/func.png')

    plt.show()