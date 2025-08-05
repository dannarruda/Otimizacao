#Nesse código, foi necessário iniciar uma Otimização de uma planta integrada para produção de etanol a partir de uma matéria prima.
#O programa tem como função maximizar o lucro operacional a partir de variáveis de decisão, tais como: Temperatura de reação e tempo de residência no reator, temperatura
#do trocador e refluxo da coluna
#A Otimização contou com restrições das variáveis
#Foi usada uma subrotina - scipy.optimize


from numpy import exp
from scipy.optimize import minimize


#Função Producao de Etanol
def producao_de_etanol(x):
    Tr = x[0]                       #Temperatura de reação
    tr = x[1]                      #Tempo de residência
    PetOH = (0.45 * (1 - exp(-0.1 * (Tr - 300))) * (1 - exp(-0.5 * tr)) * 1000)
    return PetOH

#Função Receita
def receita(PetOH):
    R = 2.0 * PetOH
    return R

#Função Custo de Matérias Primas
def custo_materia_prima():
    return 0.5*(1.2)*1000

#Função Custo de Energia
def custo_de_energia(x):
    Tth = 366                         #Temperatura de saída do trocador de calor
    Re = 1                            # Refluxo na coluna de destilação
    CE = 0.02*((Tth-(300+0.325))**2+0.1*Re**2)
    return CE

#Função Penalidade
def penalidade(x):
    Re = x[3]
    return 0.05 * Re * 1000

# Função Objetivo
def objetivo(x):
    PetOH = producao_de_etanol(x)
    R = receita(PetOH)
    CMP = custo_materia_prima()
    CE = custo_de_energia(x)
    P = penalidade(x)
    Lucro = R - CMP - CE - P
    return -Lucro


#Definindo as restrições de desigualdade

def restricao_Tr_min(x):
    return x[0] - 300  #Tr >= 300

def restricao_Tr_max(x):
    return 350 - x[0]  # Tr <= 350

def restricao_tr_min(x):
    return x[1] - 1   # tr >= 1

def restricao_tr_max(x):
    return 5 - x[1]    # tr <= 5

def restricao_Tth_min(x):
    return x[2] - 310  #Tr >= 310

def restricao_Tth_max(x):
    return 400 - x[2]  # Tr <= 400

def restricao_Re_min(x):
    return x[3] - 1 # Re >= 1

def restricao_Re_max(x):
    return 10 - x[3]    # tr <= 10

def restricao_Tr_Tth(x):
    return x[2] - (x[0] + 16)

def restricao_PEtOH(x):
    return producao_de_etanol(x) - 500  # PetOH >= 500

def restricao_custo_de_energia_max(x):
    return 100 - custo_de_energia(x)   # CE =< 100


# Lista de restrições
restricoes = [
    {'type': 'ineq', 'fun': restricao_Tr_min},
    {'type': 'ineq', 'fun': restricao_Tr_max},
    {'type': 'ineq', 'fun': restricao_tr_min},
    {'type': 'ineq', 'fun': restricao_tr_max},
    {'type': 'ineq', 'fun': restricao_Tth_min},
    {'type': 'ineq', 'fun': restricao_Tth_max},
    {'type': 'ineq', 'fun': restricao_Re_min},
    {'type': 'ineq', 'fun': restricao_Re_max},
    {'type': 'ineq', 'fun': restricao_Tr_Tth},
    {'type': 'ineq', 'fun': restricao_PEtOH},
    {'type': 'ineq', 'fun': restricao_custo_de_energia_max},
]


# Ponto inicial
x0 = [320, 2, 330, 5]

# Otimização
result = minimize(objetivo, x0, method='SLSQP', constraints=restricoes)

# Resultados
print("Sucesso:", result.success)
print("x ótimo:", result.x)
print("Lucro Máximo:", -result.fun)
print("Mensagem do otimizador:", result.message)
