##################################################################|TRABALHO FINAL - OTIMIZAÇÃO NÃO LINEAR|#########################################################################################
# Universidade Estadual de Montes claros - UNIMONTES
# Curso de Engenharia de Sistema - Otimização Não Linear - 5° Período
# Professor - Lenir De Abreu Júnior
#
# Equipe 01 — Métodos de Direção de Busca - Direções Aleatória
#   == GRUPO ==
#   -> Mateus Moreira
#   -> Uedson Gaiek
#   -> Maria Eduarda
#   -> Thiago
###################################################################################################################################################################################################
import sympy as sym
from os import system
from random import uniform

# == Classe para definir a condição de parada ==
class CondParada:
    def __init__(self, lmt) -> None:
        self.limite = lmt
        self.parada = False
    
    # 01 Condição de parada baseada no Grandiente. O gradiente sempre é zero no ponto de máximo
    def conferir_parada_grad(self, grd):
        if (abs(grd[0]) < self.limite or abs(grd[1]) < self.limite) : self.parada = True

    # 02

    # 03

# == Função principal ==
def otimz_rand(interval):
    # Definições de variáveis
    alpha, x1, x2 = sym.symbols('a x₁ x₂') # xi = xᵢ
    v = [4, 4] # X₀
    f = (((x1-3)**2)/4) + (((x2-2)**2)/9) + 13 # Função
    k = 0 # N° de iterações

    # Condição de parada
    cond = CondParada(0.001) # Tolerância de 0.001 (0.1%)

    # Loop principal
    while (not cond.parada) and (k < 300):
        d = [uniform(interval[0], interval[1]), uniform(interval[0], interval[1])] # Pontos aleatórios entre n e m
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

        # Confere a condição de parada
        cond.conferir_parada_grad(d)
        k += 1
    
    # Resultado final
    print(f'Ponto ótimo X* = ({v[0]:.4f}, {v[1]:.4f}) | {k} iterações')

    # Plot do gráfico

# Roda o programa
if __name__ == '__main__':
    system('cls')
    otimz_rand([-1, 1])