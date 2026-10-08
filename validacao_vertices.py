# -*- coding: utf-8 -*-
"""
Validação dos resultados do artigo (seções 4.2, 4.3 e 4.6; Tabelas 5, 6, 12 e 13).

O método diferencial aproxima f em torno de r por um polinômio de Taylor. Este script mede
quanto essa aproximação erra, comparando as estimativas de erro máximo com o máximo EXATO,
obtido avaliando f nos 2^n vértices do hipercubo r ± Δr.

  Parte 1 (seção 4.2): derivadas simbólicas do SymPy x formas fechadas (50 dígitos) e
                       x diferenças finitas centrais  (Modelo Exemplo, gás ideal).
  Parte 2 (4.3 e 4.6): para P2, W (gás ideal), ET0 (FAO-56) e q (parede de 3 camadas):
                       Me (1ª ordem), Me2 (2ª ordem) e máximo exato, para Δp = 1, 2, 5 e 10 %,
                       além da variação individual de +5 % em cada variável.

Uso:      python validacao_vertices.py
Requer:   sympy (instala o mpmath junto).      Tempo: cerca de 10 a 30 segundos.
Os valores nominais dos modelos FAO-56 e parede são ilustrativos (não são medições).
"""
import itertools
import sympy as sym
import mpmath as mp


def fmt(v, d=5):
    return f"{v:.{d}f}".replace(".", ",")


# ------------------------------------------------------------------ modelos
def modelo_gas():
    x = sym.symbols("P1 V1 V2 cp cv", real=True)
    P1, V1, V2, cp, cv = x
    g = cp / cv
    P2 = P1 * (V1 / V2) ** g
    return x, {"P2": P2, "W": (P1 * V1 - P2 * V2) / (g - 1)}, [101.0, 0.1, 1.6, 1040.0, 740.0]


def modelo_fao56():
    x = sym.symbols("Tmax Tmin RHmax RHmin u2 Rn G P", real=True)
    Tmax, Tmin, RHmax, RHmin, u2, Rn, G, P = x
    e0 = lambda T: 0.6108 * sym.exp(17.27 * T / (T + 237.3))
    Tm = (Tmax + Tmin) / 2
    es = (e0(Tmax) + e0(Tmin)) / 2
    ea = (e0(Tmin) * RHmax / 100 + e0(Tmax) * RHmin / 100) / 2
    D = 4098 * e0(Tm) / (Tm + 237.3) ** 2
    g = 0.665e-3 * P
    ET0 = (0.408 * D * (Rn - G) + g * 900 / (Tm + 273) * u2 * (es - ea)) / (D + g * (1 + 0.34 * u2))
    return x, {"ET0": ET0}, [30.0, 18.0, 85.0, 45.0, 2.0, 15.0, 0.6, 95.0]


def modelo_parede():
    x = sym.symbols("Tin Tout hin hout L1 L2 L3 k1 k2 k3", real=True)
    Tin, Tout, hin, hout, L1, L2, L3, k1, k2, k3 = x
    R = 1 / hin + L1 / k1 + L2 / k2 + L3 / k3 + 1 / hout
    return x, {"q": (Tin - Tout) / R}, [40.0, 10.0, 8.0, 25.0, 0.02, 0.15, 0.05, 0.7, 1.0, 0.04]


# ------------------------------------------------------------------ parte 1: seção 4.2
def parte1():
    print("=" * 78)
    print("PARTE 1 - Validação das derivadas (seção 4.2), Modelo Exemplo P2 = P1 (V1/V2)^(cp/cv)")
    x, saidas, r = modelo_gas()
    P2 = saidas["P2"]
    val = dict(zip(x, r))
    As = [float(sym.diff(P2, xi).subs(val)) for xi in x]          # derivadas simbólicas, em precisão dupla

    mp.mp.dps = 50                                                  # formas fechadas, 50 dígitos
    R = [mp.mpf(repr(v)) for v in r]
    g = R[3] / R[4]; q = R[1] / R[2]; fm = R[0] * q ** g; L = mp.log(q)
    fechadas = [q ** g, g * fm / R[1], -g * fm / R[2], fm * L / R[4], -fm * L * R[3] / R[4] ** 2]
    dmax = max(abs(mp.mpf(a) - b) / abs(b) for a, b in zip(As, fechadas))
    print(f"  Simbólico (double) x formas fechadas (50 dígitos): desvio relativo máximo = {float(dmax):.2e}")

    F = lambda v: v[0] * (v[1] / v[2]) ** (v[3] / v[4])
    print("  Diferenças finitas centrais, passo h = eta * r:")
    for eta in (1e-3, 1e-4, 1e-5, 1e-6):
        e = 0.0
        for j in range(5):
            h = eta * r[j]; a = list(r); b = list(r); a[j] += h; b[j] -= h
            fd = (F(a) - F(b)) / (2 * h)
            e = max(e, abs(fd - float(fechadas[j])) / abs(float(fechadas[j])))
        print(f"     eta = {eta:g}:  desvio relativo máximo = {e:.2e}")


# ------------------------------------------------------------------ parte 2: seções 4.3 e 4.6
def parte2():
    print("\n" + "=" * 78)
    print("PARTE 2 - Primeira e segunda ordem x máximo exato nos vértices de r ± Δr")
    for titulo, (x, saidas, r) in (("Gás ideal (n = 5)", modelo_gas()), ("FAO-56 (n = 8)", modelo_fao56()),
                                   ("Parede de 3 camadas (n = 10)", modelo_parede())):
        n = len(r); val = dict(zip(x, r))
        for nome, F in saidas.items():
            f = sym.lambdify(x, F, "math")
            f0 = f(*r)
            As = [sym.lambdify(x, sym.diff(F, xi), "math")(*r) for xi in x]
            H = [[0.0] * n for _ in range(n)]
            for i in range(n):
                for j in range(i, n):
                    H[i][j] = H[j][i] = sym.lambdify(x, sym.diff(F, x[i], x[j]), "math")(*r)
            print(f"\n--- {titulo}: saída {nome},  f(r) = {fmt(f0)} ---")
            print("   Δp     Me (1ª)   Me2 (2ª)   máx. exato   Me/exato   Me2/exato   (Tabelas 6 e 12)")
            for dp in (0.01, 0.02, 0.05, 0.10):
                dr = [dp * v for v in r]
                Me1 = sum(abs(a) * abs(d) for a, d in zip(As, dr))
                Me2 = Me1 + 0.5 * sum(abs(H[i][j]) * abs(dr[i]) * abs(dr[j]) for i in range(n) for j in range(n))
                ex = max(abs(f(*[r[i] + s[i] * dr[i] for i in range(n)]) - f0)
                         for s in itertools.product((-1, 1), repeat=n))
                print(f"  {dp*100:4.0f} %  {fmt(Me1):>9} {fmt(Me2):>10} {fmt(ex):>12}  {fmt(100*Me1/ex,1):>8} %  {fmt(100*Me2/ex,1):>8} %")
            if titulo.startswith("Gás"):
                dr = [0.05 * v for v in r]
                print("  Variação individual de +5 % em cada variável (Tabelas 5 e 13):")
                print("     variável      exata     1ª ordem     2ª ordem   desvio 1ª   desvio 2ª")
                for v in range(n):
                    t = list(r); t[v] += dr[v]; ex = f(*t) - f0
                    p1 = As[v] * dr[v]; p2 = p1 + 0.5 * H[v][v] * dr[v] ** 2
                    d1 = abs(p1 - ex) / abs(ex) * 100; d2 = abs(p2 - ex) / abs(ex) * 100
                    print(f"     {str(x[v]):<8} {fmt(ex):>10} {fmt(p1):>12} {fmt(p2):>12} {fmt(d1,1):>9} % {fmt(d2,1):>9} %")
            # escalonamento uniforme de +5 %
            unif = f(*[1.05 * v for v in r]) / f0 - 1
            print(f"  Escalonamento uniforme de +5 % em todas as variáveis: variação relativa exata = {fmt(unif)}")


if __name__ == "__main__":
    parte1()
    parte2()
    print("\nConcluído.")
