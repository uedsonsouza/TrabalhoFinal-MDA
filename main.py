##################################################################|TRABALHO FINAL - OTIMIZAÇÃO NÃO LINEAR|#########################################################################################
# Universidade Estadual de Montes Claros - UNIMONTES
# Curso de Engenharia de Sistema - Otimização Não Linear - 5° Período
# Professor - Lenir De Abreu Júnior
#
# Equipe 01 — Métodos de Direção de Busca - Direções Aleatória
#   == GRUPO ==
#   -> Mateus Moreira Durães
#   -> Uedson Gaiek
#   -> Maria Eduarda Alves Cardoso
#   -> Thiago Emanoel Brito Santos
###################################################################################################################################################################################################
import sympy as sym
import matplotlib.pyplot as plt
from os import system
from random import uniform

# == Configurações ==
STOP_MODE = 1 # Modo de condição de parada [0: Gradiente, 1: Variação, 2: ??]
MAX = 300     # Número máximo de iterações

# Retorna um vetor gradiente de uma função
def grad(f, arg1, arg2, ptn):
    a = sym.derivative_f = f.diff(arg1) # Derivada parcial em função de x₁
    b = sym.derivative_f = f.diff(arg2) # Derivada parcial em função de x₂

    arr = [ # Vetor Gradiente
    -( (a.subs(arg1, ptn[0])).subs(arg2, ptn[1]) ),
    -( (b.subs(arg1, ptn[0])).subs(arg2, ptn[1]) )
    ]

    return arr

# == Classe para definir a condição de parada ==
class CondParada:
    def __init__(self, lmt) -> None:
        self.limite = lmt
        self.mode = STOP_MODE
        self.parada = False
    
    # TODO >> Conferir últimos 5 pontos ao invés de apenas o último 
    # 01 Condição de parada baseada no Grandiente. O gradiente sempre é zero no ponto de máximo
    def conferir_parada_grad(self, grd):
        if (abs(grd[0]) < self.limite or abs(grd[1]) < self.limite) : return True

    # 02 Condição de parada baseada na Variação entre os últimos 5 pontos
    def conferir_parada_var(self, arg):
        if len(arg) < 5 : return False
        for i, j in zip(range(-5, -1), range(-4, -1)):
            if (abs(arg[i][0] - arg[j][0]) > self.limite) or (abs(arg[i][1] - arg[j][1]) > self.limite): return False
        return True

    # 03

    #
    def set_stop(self, grd = [], arg = []):
        if   (self.mode == 0) and (self.conferir_parada_grad(grd)): self.parada = True
        elif (self.mode == 1) and (self.conferir_parada_var(arg)): self.parada = True
       #elif (self.mode == 2) and 

# == Função principal ==
def otimz_rand(interval):
    # Definições de variáveis
    alpha, x1, x2 = sym.symbols('a x₁ x₂') # xi = xᵢ
    v = [4, 4] # X₀
    f = (((x1-3)**2)/4) + (((x2-2)**2)/9) + 13 # Função
    k = 0 # N° de iterações

    # Armazenando informações
    pnts = [v]

    # Condição de parada
    cond = CondParada(0.001) # Tolerância de 0.1%

    print('Calculando. Por favor espere...')
    # Loop principal
    while (not cond.parada) and k < MAX:
        # TODO >> Mudar para variáveis aleatórias nomais ao invés de uniformes
        d = [uniform(interval[0], interval[1]), uniform(interval[0], interval[1])] # Pontos aleatórios entre n e m OBS: Mudar para normal

        p = [ # Pontos com relação de um alpha (direção)
            (v[0] + alpha * d[0]),
            (v[1] + alpha * d[1])
        ]

        temp = (f.subs(x1, p[0])).subs(x2, p[1]) # Substitui os pontos na função
        f_dif = temp.diff(alpha) # Diferencia a nova função em relação a alpha
        res = sym.solve(f_dif, alpha) # Resolve a função para f_dif = 0 para achar o valor de alpha

        if len(res) > 0:
            v = [ # Novos pontos são atribuídos
                p[0].subs(alpha, res[0]),
                p[1].subs(alpha, res[0])
            ]
        else: v = p
        pnts.append(v)

        # Confere a condição de parada
        cond.set_stop( grd = grad(f, x1, x2, v), arg = pnts)
        k += 1

        # TODO >> Atribuir o break no loop while
        # Abortar após grande número de iterações
        #if k >= MAX: break
    
    # Resultado final
    system('cls')
    print(f'{pnts}')
    print(f'Ponto ótimo X* = ({v[0]}, {v[1]}) | {k} iterações')

    # TODO >> Plot em 3D
    # Plot do gráfico
    x = []
    for x_i in range(len(pnts)): x.append(x_i)

    figure, axis = plt.subplots(2, 1)
    axis[0].plot(x, pnts)
    axis[0].set_title("f(x₁, x₂) x Iteração")

    plt.show()
    



# Roda o programa
if __name__ == '__main__':
    system('cls')
    otimz_rand([-1, 1])