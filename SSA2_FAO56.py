# Modelo de Penman-Monteith (FAO-56), n = 8 variáveis (Tabela 10 do artigo)
# Requer Relatorio2.py na mesma pasta.  Execute:  python SSA2_FAO56.py
# Gera ssa2_fao56.html  (valores nominais ilustrativos, Δp = 5 % em todas as variáveis)

import sympy as sym
from   Relatorio2 import relatorio2

n = 8
r = [0.0]*n;  Δp = [0.0]*n

r[0] =    30.0;   Δp[0] = 0.05   # Tmax (°C)
r[1] =    18.0;   Δp[1] = 0.05   # Tmin (°C)
r[2] =    85.0;   Δp[2] = 0.05   # RHmax (%)
r[3] =    45.0;   Δp[3] = 0.05   # RHmin (%)
r[4] =     2.0;   Δp[4] = 0.05   # u2 (m/s)
r[5] =    15.0;   Δp[5] = 0.05   # Rn (MJ m-2 dia-1)
r[6] =     0.6;   Δp[6] = 0.05   # G (MJ m-2 dia-1)
r[7] =    95.0;   Δp[7] = 0.05   # P (kPa)

var = ['Tmax','Tmin','RHmax','RHmin','u2','Rn','G','P']
x = list(sym.symbols(var, real=True))
Tmax, Tmin, RHmax, RHmin, u2, Rn, G, P = x

e0 = lambda T: 0.6108*sym.exp(17.27*T/(T + 237.3))
Tm = (Tmax + Tmin)/2
es = (e0(Tmax) + e0(Tmin))/2
ea = (e0(Tmin)*RHmax/100 + e0(Tmax)*RHmin/100)/2
Δ  = 4098*e0(Tm)/(Tm + 237.3)**2
γ  = 0.665e-3*P
F  = sym.Matrix([(0.408*Δ*(Rn - G) + γ*900/(Tm + 273)*u2*(es - ea))/(Δ + γ*(1 + 0.34*u2))])
nomes = ['ET0']

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

relatorio2("ssa2_fao56.html", var, r, Δp, Δr, saidas)
