from sympy import *
import numpy as np
from scipy.optimize import fmin

x1, x2, a = symbols("x1 x2 a")
X = [(3, 4),]
rosenbrock = 100*(x2-(x1**2))**2 + (x1-1)**2
restrictions = (x1/9)+(x2/16)
fx= rosenbrock+restrictions
print (fx)

x1/9 + x2/16 + (x1 - 1)**2 + 100*(-x1**2 + x2)**2

def gradiente(fx, xk):
    g = []
    g.append(diff(fx, x1))
    g.append(diff(fx, x2))

    for i in range(len(g)):
        g[i] = g[i].subs(([(x1, xk[0]), (x2, xk[1])]))

    return g

def f(m):
    return X[k] + a *D[k]

k=0
G = []
D = []
A = []

while k<1:
    print("-"*30, f"K = {k}", "-"*30)
    G.append(gradiente(fx, X[k]))
    print(f"\nGradiente: {G[k]}")

    dk = [G[k][0] * -1, G[k][1] * -1]
    D.append(dk)
    print(f"\nDireção dk: {D[k]}")

    Aux = [X[k][0] + a * dk[0], X[k][1] + a * dk[1]]
    print(f"\nValor da função f(xk+adk): {Aux}")

    fa = fx.subs([(x1, Aux[0]), (x2, Aux[1])], simplify = False)
    fa = simplify(fa)
    print(f"\nF(a): {fa}")

    dadf = diff(fa, a)
    print(f"\nDerivada F(a): {dadf}")

    solution_dadf = fmin(f, 0, disp=0)
    print(f"\nSolucão Da/Df: {solution_dadf}")

    k+=1
    print("-"*70)
