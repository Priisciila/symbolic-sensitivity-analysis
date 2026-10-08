import sympy as sym
from   Relatorio import relatorio

n = 5
r = [0.0]*n;  Δp = [0.0]*n

r[0] =  101.0;   Δp[0] = 0.05
r[1] =    0.1;   Δp[1] = 0.05
r[2] =    1.6;   Δp[2] = 0.05
r[3] = 1040.0;   Δp[3] = 0.05
r[4] =  740.0;   Δp[4] = 0.05

var = ['P1','V1','V2','cp','cv']
P1, V1, V2, cp, cv = sym.symbols( " ".join(str(e) for e in var), real=True )
x = [sym.var(v) for v in var]

f = P1*(V1/V2)**(cp/cv)

Δr = [0.0]*n; S  = [0.0]*5
As = [0.0]*n; Rs = [0.0]*n; Pie = [0.0]*n; Me = [0.0]*n; Re = [0.0]*n; Mre = [0.0]*n

val = {}
for v,i in zip(x,r): val[v] = i

for v in range(0,n):
    Δr [v] = Δp[v]*r[v]
    As [v] = sym.diff(f,x[v])
    Rs [v] = As[v] * r[v]  / f                ; S[0] += Rs [v]
    Pie[v] = As[v] * Δr[v]                    ; S[1] += Pie[v]
    Me [v] = sym.Abs(As[v]) * sym.Abs(Δr[v])  ; S[2] += Me [v]
    Re [v] = As[v] * Δr[v] / f                ; S[3] += Re [v]
    Mre[v] = sym.Abs(As[v]) * sym.Abs(Δr[v]/f); S[4] += Mre[v]

relatorio("ssa.html", x, r, Δp, Δr, As, Rs, Pie, Me, Re, Mre, S, val)
