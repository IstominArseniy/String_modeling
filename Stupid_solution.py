from matplotlib import pyplot as plt
import math
import scipy.optimize
import numpy as np

# fitting numbers
M1 = 4
M2 = 3
# all in SI
rho_0 = 1  # density of air
mu = 1.45e-5  # dynamic viscosity of air
Qt = 18000  # specific constant in thermoelasticity
#  string params
rho = 2500  # string dencity
r = 0.44e-3  # radius of string
T = 102.6  # Tension
E = 8.6e9  # Young module
S = math.pi * r ** 2  # cross area
c = (T / rho / S) ** 0.5  # wave speed
eta = 9.5e-4
I = math.pi * r ** 4 / 4


def beta(w):
    """
    dispersion function
    :param w: circular frequency
    :return: wave vector
    """
    k2 = E * I / S / rho
    return (((c ** 4 + 4 * w ** 2 * k2) ** 0.5 - c ** 2) / 2 / k2) ** 0.5


def Za(f):
    w = math.pi * 2 * f
    return 2 * math.pi * (mu + r * (w * 2 * mu * rho_0) ** 0.5)


def Zt(f):
    w = math.pi * 2 * f
    return rho * math.pi * r ** 2 * w / Qt


def Zv(f):
    w = math.pi * 2 * f
    return eta * E * I / c * (beta(w)) ** 3


def Zth(f):
    return Zv(f) + Za(f) + Zt(f)


def Zfoster(f, a11, a12, a13, a14, b11, b12, b13, b14, a21, a22, a23,  b21, b22, b23):
    a1 = [a11, a12, a13, a14]
    b1 = [b11, b12, b13, b14]
    a2 = [a21, a22, a23]
    b2 = [b21, b22, b23]
    ans = 0.0
    ans1 = 0.0
    ans2 = 0.0
    w = 2 * math.pi * f
    for ind in range(M1):
        ans1 += b1[ind] / (a1[ind] ** 2 + w ** 2)
    for ind in range(M2):
        ans2 += b2[ind] / (a2[ind] ** 2 + w ** 2)
    ans = w ** 2 * ans1 + w ** 2 * beta(w) ** 2 * ans2
    # ans = params[0] * f + params[1]
    return ans


# freq = np.arange(100, 20000, 50)
freq = np.logspace(2, 4, 500)
p0 = [0, 141.95, 725.57, 3340.6, 2.355e-4, 2.034e-4, 4.674e-4, 1.199e-3, 0, 48832, 193909, 1.11e-8, 2.13e-7, 5.78e-7]
ppot, pcov = scipy.optimize.curve_fit(Zfoster, freq, Zth(freq), maxfev=100000, p0=p0)
# print(ppot)
print(Zfoster(5000, *p0))
plt.yscale("log")
plt.xscale("log")
plt.plot(freq, Za(freq), color='blue')
plt.plot(freq, Zt(freq), color='red')
plt.plot(freq, Zv(freq), color='green')
plt.plot(freq, Zth(freq), color='black')
plt.plot(freq, Zfoster(freq, *ppot), color='yellow')
plt.show()
