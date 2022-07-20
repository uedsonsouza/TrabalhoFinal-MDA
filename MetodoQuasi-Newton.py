from builtins import print
from sympy import *
import numpy as np
import matplotlib.pyplot as plt

#função para calcular o gradiente de uma função de variáveis simbolicas
def gradiente(fx, symx):
    g = []
    for x in symx:
        g.append(diff(fx, x))
    return(g)

#função para substituir o ponto x0 = (x1, x2, ...) nas coordenadas do gradiente
def subsGrad(gx, symx, x0):
    subsFx = gx[:]
    for i in range(len(subsFx)):
        for j in range(len(symx)):
            subsFx[i] = subsFx[i].subs(symx[j], x0[j])
        
    return subsFx

def calculaCDFP(vk_temp, rk_temp, hk):
    # ((vk*vk')/(vk'*rk))-((hk*rk*rk'*hk)/(rk'*hk*rk))
    
    #transformando em vetores numpy
    vk = np.array(vk_temp)
    rk = np.array(rk_temp)

    #numerador 1ª parte
    parte1 = np.outer(vk, vk.T)
    #denominador 1ª parte
    parte2 = np.dot(vk.T, rk)
    
    #divisão 1ª parte
    parte3 = parte1 / parte2

    #numerador 2ª parte
    parte4 = np.dot(hk, rk)
    parte4 = np.reshape(parte4, (len(rk), 1))
    parte5 = np.dot(rk.T, hk)
    parte6 = np.outer(parte4, parte5)

    #denominador 2ª parte
    parte7 = np.dot(rk.T, hk)
    parte8 = np.dot(parte7, rk)
    

    #divisão 2ª parte
    parte9 = parte6 / parte8

    #subtração dos dois membros
    CDFP = parte3 - parte9

    return CDFP

def calculaCBFGS(vk_temp, rk_temp, hk):
    # (1+((rk'*hk*rk)/(rk'*vk)))*((vk*vk')/(vk'*rk))-((vk*rk'*hk + hk*rk*vk')/(rk'*vk))

    #transformando em vetores numpy
    vk = np.array(vk_temp)
    rk = np.array(rk_temp)

    #numerador primeira parte
    parte1 = np.dot(rk.T, hk)
    parte2 = np.dot(parte1, rk)
    #denomiador primeira parte
    parte3 = np.dot(rk.T, vk)

    #divisao primeira parte
    parte4 = parte2 / parte3
    
    #soma primeira parte
    parte5 = 1 + parte4

    #numerador segunda parte
    parte6 = np.outer(vk, vk.T)
    #denominador segunda parte
    parte7 = np.dot(vk.T, rk)
    
    #divisao segunda parte
    parte8 = parte6 / parte7

    #primeira * segunda parte
    parte9 = parte5 * parte8

    #numerador terceira parte
    parte10 = np.outer(vk, rk.T)
    parte11 = np.dot(parte10, hk)

    parte12 = np.dot(hk, rk)
    parte12 = np.reshape(parte12, (len(rk), 1))
    parte13 = parte12 * vk.T

    parte14 = parte11 + parte13

    #denominador terceira parte
    parte15 = np.dot(rk.T, vk)

    #divisao terceira parte
    parte16 = parte14 / parte15

    BFGS = parte9 - parte16
    
    return BFGS

def critParada(gk, funcoes):
    if len(funcoes) < 6:
        elevQuadrado = list(map(lambda x: x ** 2, gk))
        soma = sum(elevQuadrado)
        raiz = sqrt(soma)

        return raiz > 0.001
    else:
        try:
            Deltaf = max(funcoes) - min(funcoes)
            fcincomais = max(funcoes[-6:])
            fcincomenos = min(funcoes[-6:])
            deltinhaf = fcincomais - fcincomenos
            if deltinhaf < (0.0001 * Deltaf):
                return False
            else:
                return True

        except:
            elevQuadrado = list(map(lambda x: x ** 2, gk))
            soma = sum(elevQuadrado)
            raiz = sqrt(soma)

            return raiz > 0.001

def subsFx(fx, symx, x):
    subsFx = fx
    for i in range(len(symx)):
        subsFx = subsFx.subs(symx[i], x[i])

    return subsFx

def QuasiNewton(fx, x, symx, n):

    fg = gradiente(fx, symx)
    hk = np.eye(n)
    gk = subsGrad(fg,symx, x[0])
    k = 0
    funcoes = []
    a = symbols("a")
    while critParada(gk, funcoes):
        X = x[k] - a * np.dot(hk, gk)
        
        Fdealfa = fx.subs([(symx[0], X[0]), (symx[1], X[1])])

        DifFdealfa = diff(Fdealfa, a)

        try:
            solveAlfa = nroots(DifFdealfa, maxsteps = 10000)
        except:
            print("\n\nNão foi possível encontrar o ponto ótimo desta função")
            exit(1)

        alfa = solveAlfa[0]
        
        try:
            for i in range(len(solveAlfa)):
                if(Fdealfa.subs(a, solveAlfa[i]) < Fdealfa.subs(a, alfa)):
                    alfa = solveAlfa[i]
        except:
            alfa = solveAlfa[0]

        for i in range(len(X)):
            X[i] = X[i].subs(a, alfa)


        G = subsGrad(fg, symx, X)

        vk = list(map(lambda i, j: i - j, x[k], X))
        
        rk = list(map(lambda i, j: i - j, gk, G))

        CDFP = calculaCDFP(vk, rk, hk)
        CBFGS = calculaCBFGS(vk, rk, hk)

        beta = np.random.rand()

        ck = (1-beta) * CDFP + beta * CBFGS

        H = hk + ck;

        k += 1

        x.append(X)

        gk = np.copy(G)
        
        funcoes.append(subsFx(fx, symx, X))

        hk = np.copy(H)

    print(f"\n\nQuantidade de execuções: {len(x)}")
    print(f"Ponto ótimo encontrado: {tuple(x[-1])}")
    print(f"A função no ponto: {fx.subs([(symx[0], x[-1][0]), (symx[1], x[-1][1])])}")

def graph3D(fx, x, symx):
    plt.rcParams['figure.figsize'] = 8, 6
    plotting.plot3d(fx, title = f"Gráfico da Função: {fx}")
    
def graphNivel(fx, x, symx):
    fig, ax = plt.subplots(figsize=(102, 10))
    #curvas de nível
    xvec = np.linspace(-5, int(max(x[:][0])) + 5, 500)
    yvec = np.linspace(-5, int(max(x[:][1])) + 5, 500)
    xgraf, ygraf = np.meshgrid(xvec, yvec)

    funclam = lambdify(symx, fx)
    funcao = funclam(xgraf, ygraf)

    contornp = ax.contour(xgraf, ygraf, funcao, 100)
    plt.grid()


    plt.scatter(x[-1][0], x[-1][1], s=100, c="red")
    plt.scatter(x[0][0], x[0][1], s=100, c="green")


    for i in range(len(x)-1):
        plt.plot([x[i][0], x[i+1][0]], [x[i][1], x[i+1][1]], linewidth = 3, color = "green")
        plt.scatter(x[i][0], x[i][1], s=50, c="green")


    plt.annotate(
            f"x0: {tuple(x[0])}",
            xy=(x[0][0], x[0][1]), xytext=(-20, 20),
            textcoords='offset points', ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.5', fc='yellow'),
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'), fontsize=13)

    plt.annotate(
            f"x*: ({x[-1][0]:.2f}, {x[-1][1]:.2f})",
            xy=(x[-1][0], x[-1][1]), xytext=(-20, 20),
            textcoords='offset points', ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.5', fc='yellow'),
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'), fontsize=13)

    plt.title("Gráfico de busca - Quasi Newton")
    plt.show()

def graphConvertion(fx, x, symx):
    iteracoes = []
    valores = []

    for i in range(0, len(x)):
        iteracoes.append(i)
        valores.append(fx.subs([(symx[0], x[i][0]), (symx[1], x[i][1])]))

    plt.figure(figsize=(10, 8))
    plt.xlabel("Iteracoes", fontsize = 16)
    plt.ylabel("F(x)", fontsize = 16)
    plt.title("Gráfico de convergência", fontsize = 16)
    plt.plot(iteracoes, valores)
    plt.yscale("log")
    plt.xticks(iteracoes)
    plt.show()

if __name__ == '__main__':
    n = int(input("Digite quantas variáveis sua função possui: "))
    symx = symbols(f"x(1:{n+1})")

    f = str(input("Digite a sua função: "))
    fx = sympify(f, convert_xor=True)

    x0 = list(map(float,input("Digite o seu ponto x0 (ex: 1 2): ").strip().split()))[:n]

    x = [x0, ]

    QuasiNewton(fx, x, symx, n)


    resp = str(input("\n\nGostaria de visualizar os gráficos?(sim/nao) "))


    if resp in ["sim", "Sim", "SIM", "s", "S"]:
        graph3D(fx, x, symx)

        graphNivel(fx, x, symx)

        graphConvertion(fx, x, symx)


