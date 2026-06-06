import numpy as np
import matplotlib.pyplot as plt

def numerical_diff(f,x):
    h = 1e-4 # 0.0001
    nd = (f(x+h) - f(x-h))/(2 * h)
    return nd

def func_1(x):
    return x**2

x = np.arange(0.0,2.0,0.01) # 0から20まで0.1刻みのベクトルを生成

y = func_1(x)

y_diff = [numerical_diff(func_1, xi) for xi in x]

plt.xlabel("x")
plt.ylabel("f(x)")
plt.plot(x,y)
plt.plot(x,y_diff)

plt.show()