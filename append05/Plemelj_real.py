import numpy as np
import matplotlib.pyplot as plt

def func(a, x_list):
    y_list = [a/(a**2 + x**2) for x in x_list]
    return y_list

if __name__=='__main__':
    a0 = 1
    a1 = 0.1
    a2 = 0.01
    param0 = 'alpha = {0}'.format(a0)
    param1 = 'alpha = {0}'.format(a1)
    param2 = 'alpha = {0}'.format(a2)
    x = np.linspace(-3, 3, 1000)
    y0 = func(a0, x)
    y1 = func(a1, x)
    y2 = func(a2, x)

    # drawin graphs
    fig = plt.figure(tight_layout=True)
    ax = fig.add_subplot(111)
    ax.grid()
    ax.plot(x, y0, label=param0)
    ax.plot(x, y1, label=param1)
    ax.plot(x, y2, label=param2)
    ax.set_xlabel('omega')
    ax.set_ylim(0,20)
    plt.legend(loc='upper left')

    fig.savefig('./png/func.png')

    plt.show()