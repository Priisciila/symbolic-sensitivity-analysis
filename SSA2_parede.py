# Parede composta de 3 camadas, n = 10 variáveis (Tabela 11 do artigo)
# Requer Relatorio2.py na mesma pasta.  Execute:  python SSA2_parede.py
# Gera ssa2_parede.html  (valores nominais ilustrativos, Δp = 5 % em todas as variáveis)

import sympy as sym
from   Relatorio2 import relatorio2

n = 10
r = [0.0]*n;  Δp = [0.0]*n

r[0] =   40.00;   Δp[0] = 0.05   # Tin (°C)
r[1] =   10.00;   Δp[1] = 0.05   # Tout (°C)
r[2] =    8.00;   Δp[2] = 0.05   # hin (W m-2 K-1)
r[3] =   25.00;   Δp[3] = 0.05   # hout (W m-2 K-1)
r[4] =    0.02;   Δp[4] = 0.05   # L1 reboco (m)
r[5] =    0.15;   Δp[5] = 0.05   # L2 concreto (m)
r[6] =    0.05;   Δp[6] = 0.05   # L3 EPS (m)
r[7] =    0.70;   Δp[7] = 0.05   # k1 (W m-1 K-1)
r[8] =    1.00;   Δp[8] = 0.05   # k2 (W m-1 K-1)
r[9] =    0.04;   Δp[9] = 0.05   # k3 (W m-1 K-1)

var = ['Tin','Tout','hin','hout','L1','L2','L3','k1','k2','k3']
x = list(sym.symbols(var, real=True))
Tin, Tout, hin, hout, L1, L2, L3, k1, k2, k3 = x

R  = 1/hin + L1/k1 + L2/k2 + L3/k3 + 1/hout
F  = sym.Matrix([(Tin - Tout)/R])
nomes = ['q']

val = dict(zip(x,r))
Δr  = [Δp[v]*r[v] for v in range(n)]
J   = F.jacobian(x)
num = lambda e: float(e.subs(val))
saidas = []

for i in range(len(F)):
    f0 = num(F[i]);  H = sym.hessian(F[i], x)
    As = [num(J[i,v]) for v in range(n)]
    Hs = [[num(H[v,w]) for w in range(n)] for v in range(n)]
    q2 = sum(Hs[v][w]*Δr[v]*Δr[w]      for v in range(n) for w in range(n))/2
    a2 = sum(abs(Hs[v][w]*Δr[v]*Δr[w]) for v in range(n) for w in range(n))/2
    Rs  = [As[v]*r[v]/f0 for v in range(n)]
    Pie = [As[v]*Δr[v]   for v in range(n)]
    Me  = [abs(As[v]*Δr[v]) for v in range(n)]
    Re  = [Pie[v]/f0 for v in range(n)]
    Mre = [Me[v]/f0  for v in range(n)]
    saidas.append(dict(nome=nomes[i], f0=f0, As=As, Rs=Rs, Pie=Pie, Me=Me,
                       Re=Re, Mre=Mre, Pie2=sum(Pie)+q2, Me2=sum(Me)+a2))

relatorio2("ssa2_parede.html", var, r, Δp, Δr, saidas)
