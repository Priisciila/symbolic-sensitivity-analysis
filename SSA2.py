import sympy as sym
from   Relatorio2 import relatorio2

n = 5
r = [0.0]*n;  Δp = [0.0]*n

r[0] =  101.0;   Δp[0] = 0.05
r[1] =    0.1;   Δp[1] = 0.05
r[2] =    1.6;   Δp[2] = 0.05
r[3] = 1040.0;   Δp[3] = 0.05
r[4] =  740.0;   Δp[4] = 0.05

var = ['P1','V1','V2','cp','cv']
x = list(sym.symbols(var, real=True))
P1, V1, V2, cp, cv = x

γ  = cp/cv
P2 = P1*(V1/V2)**γ
F  = sym.Matrix([P2, (P1*V1 - P2*V2)/(γ - 1)])
nomes = ['P2', 'W']

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

relatorio2("ssa2.html", var, r, Δp, Δr, saidas)
