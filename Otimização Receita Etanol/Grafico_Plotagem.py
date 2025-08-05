#Nesse caso, foi necessário construir o gráfico de otimização do exmeplo anterior

import numpy as np
import matplotlib.pyplot as plt
from numpy import exp

# Função Produção de Etanol
def producao_de_etanol(Tr, tr):
    return 0.45 * (1 - exp(-0.1 * (Tr - 300))) * (1 - exp(-0.5 * tr)) * 1000

# Receita
def receita(PetOH):
    return 2.0 * PetOH

# Custo de Matéria-Prima
def custo_materia_prima():
    return 0.5 * 1.2 * 1000

# Custo de Energia (Tth e Re fixos como no seu código)
def custo_de_energia():
    Tth = 366
    Re = 1
    return 0.02 * ((Tth - (300 + 0.325))**2 + 0.1 * Re**2)

# Penalidade (Re fixo)
def penalidade():
    Re = 1
    return 0.05 * Re * 1000

# Função Lucro
def lucro(Tr, tr):
    PetOH = producao_de_etanol(Tr, tr)
    R = receita(PetOH)
    CMP = custo_materia_prima()
    CE = custo_de_energia()
    P = penalidade()
    return R - CMP - CE - P

# Geração da malha de valores
Tr_vals = np.linspace(300, 350, 100)
tr_vals = np.linspace(1, 5, 100)
Tr_grid, tr_grid = np.meshgrid(Tr_vals, tr_vals)
Lucro_grid = lucro(Tr_grid, tr_grid)

# Plotagem 3D
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(Tr_grid, tr_grid, Lucro_grid, cmap='viridis', edgecolor='none')
ax.set_xlabel('Temperatura de Reação (Tr)')
ax.set_ylabel('Tempo de Residência (tr)')
ax.set_zlabel('Lucro')
ax.set_title('Lucro em função de Tr e tr')
fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)
plt.show()
