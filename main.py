import math
import numpy as np
from sympy import *
from scipy.integrate import quad
from matplotlib import pyplot as plt
from scipy.special import gamma
import pandas as pd
import math
import csv

tTheta = np.linspace(0, np.pi, 60)
theta = tTheta.reshape(60, 1)
# print(theta)
rr = np.linspace(1, 40, 39)
r = rr.reshape(1, 39)
print(r[0][1])

omega = np.cos(theta)
x = np.dot(np.sin(theta), r)
# print(omega)
R = 7.14e7
mu_0 = 4*np.pi*10**-7
init = int(1e1)
rc = 1
p = symbols("x")
q = 3e-9/(np.sqrt(2)*gamma(7.6))*pow(p, 6.5)*pow(np.e, -p/1.9)
dify = diff(q, p)
# print(dify)

def dP(t):
    return dify.subs('x', t)


def f(y):
    return dP(y)


def cni(n, i):
    return math.factorial(n) / math.factorial(i) / math.factorial(n - i)


def P_n(n,t):
    pv = 0
    # print("n is ", n)
    for m in range(0, n + 1):
        pv += cni(n + 1, m) * cni(n + 1, n - m) * (t - 1) ** (n - m) * (t + 1) ** m
        # print("m is ", m, pv)
    pv /= 2 ** n
    # print("pv is ", pv)
    return pv


def g(x, y):
    return (1 - x**2)*f(y)


# def h_n(n, x):
#     # print("h_n\n")
#     tum, ept = quad(lambda t: (1 - x ** 2) * P_n(n,t) ** 2, -1, 1)
#     # print(tum)
#     return tum


def g_n(n, x, y):
    # print("g_n\n")
    tum, ept = quad(lambda t: P_n(n, t) * g(t, y), -1, +1)
    return (n+2) / (n+1) / 8 * tum


def y_n(x, r, n, y):
    # print("y_n\n")
    res = (1 - x * x) / r
    # print(res)
    # print("test1\n")
    for v in range(0, init):
        tum_1, ept = quad(lambda t: t ** (-v - 1) * g_n(v, x, y), rc, np.inf)
        # print("test11\n")
        tum_2, ept = quad(lambda t: t ** (-v - 1) * g_n(v, x, y), rc, np.inf)
        tum_3, ept = quad(lambda t: t ** (-v - 1) * g_n(v, x, y), r, np.inf)
        res += 1 / (2 * v + 3) * (r ** (-v - 1) * (tum_1
                                                   - rc ** (2 * v + 3) * tum_2) +
                                  r ** (v + 2) * tum_3) * (1 - x ** 2) * P_n(n, x)
    # print(res)
    # print("y_n pass\n")
    return res

def writeFile(y):
    return


print("start")
y = np.random.random((60, 39)) * 100
B_r = np.zeros((60, 39))
B_theta = np.zeros((60,39))
for j in range(0, 60):
    print("j is ", j)
    for k in range(0, 39):
        print("k is ", k)
        n = 1
        tmp = y_n(x[j][k], r[0][k], n, y[j][k])
        print("y is ", y[j][k])
        print("y_ is ", tmp)
        while abs(y[j][k] - tmp) > 1e-5:
            print("n is ", n)
            if tmp < 0:
                tmp = -tmp
            y[j][k] = tmp
            n = n + 1
            tmp = y_n(x[j][k], r[0][k], n, tmp)
            print("error is ", abs(y[j][k] - tmp))
            if(n > 10):
                break
        y[j][k] = tmp

for i in range(0, 60):
    for j in range(0, 39):
        B_r[i][j] = (y[(i+1) % 60][j] - y[i][j]) / np.pi * 60 / r[1][j]**2 / np.sin(theta)[j][1]
        B_theta[i][j] = - (y[i][(j+1) % 39][j] - y[i][j]) / 1 / r[1][j] / np.sin(theta)[j][1]


print(y)



# plt.rcParams['font.sans-serif'] = ['SimHei']
#
# plt.figure(1)
# plt.xlabel("x")
# plt.ylabel("y")
# plt.semilogy(x, y)
# plt.show()




































































































































































































































































