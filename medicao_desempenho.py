# -*- coding: utf-8 -*-
"""
Comparação de métodos para obter as derivadas parciais de primeira ordem (gradiente de f em r).
Reproduz o procedimento da seção 3.5 e as Tabelas 7 a 9 do artigo.

Métodos (cada um é OPCIONAL: se o programa não estiver instalado, é pulado):
  SymPy (subs)       derivação simbólica + avaliação por substituição (como em SSA.py)
  SymPy (lambdify)   derivação simbólica + conversão em função numérica
  SymEngine          pip install symengine numpy
  Maxima             programa externo; precisa do comando "maxima" no PATH
  JAX                pip install jax              (diferenciação automática, jit, precisão dupla)
  Autograd           pip install autograd         (diferenciação automática)
  Dif. finitas       diferenças centrais, h = 1e-5 * r, com escalares NumPy (como no artigo; precisa de numpy)
  Dif. finitas (float)   o mesmo com ponto flutuante nativo do Python: 3 a 7 vezes mais rápido que a linha
                     anterior; mostra que o tempo das diferenças finitas depende da implementação

Precisão: maior desvio relativo do gradiente (precisão dupla) em relação a uma referência de
60 dígitos obtida com mpmath por diferenciação numérica (independente dos métodos testados).
Tempos (ms): "preparação" (derivação + conversão/compilação; 1ª chamada em JAX e Autograd) e
"avaliação" do gradiente completo (menor valor de várias repetições, após aquecimento).

Uso:      python medicao_desempenho.py                       (modelos padrão, ~1 a 2 min)
          python medicao_desempenho.py --modelos M1,M2,M3,W54,W104,W204
Modelos:  M1 = gás ideal (n=5), M2 = FAO-56 (n=8), M3 = parede de 3 camadas (n=10),
          W<n> = parede com K camadas, n = 4 + 2K  (ex.: W14 = 5 camadas, W204 = 100 camadas).
Opções:   --max-subs N   não executa "SymPy (subs)" para n > N (padrão 54; o custo cresce muito)
Saída:    tabelas no terminal e o arquivo resultados_desempenho.json.
Os tempos dependem do computador e das versões; use-os pela ordem de grandeza e pela ordenação.
"""
import argparse, json, math, os, platform, shutil, subprocess, tempfile, time
import sympy as sym
import mpmath as mp

mp.mp.dps = 60
try:
    import numpy as np
except Exception:
    np = None
try:
    import symengine as se
except Exception:
    se = None
try:
    import jax
    import jax.numpy as jnp
    jax.config.update("jax_enable_x64", True)
except Exception:
    jax = None
try:
    import autograd
    import autograd.numpy as anp
except Exception:
    autograd = None
MAXIMA = shutil.which("maxima")


# ------------------------------------------------------------------ modelos (uma definição serve a todos os motores)
def e0(T, xp):
    return 0.6108 * xp.exp(17.27 * T / (T + 237.3))


def M1(x, xp):
    P1, V1, V2, cp, cv = x
    return P1 * (V1 / V2) ** (cp / cv)


def M2(x, xp):                               # FAO-56 Penman-Monteith (equação diária)
    Tmax, Tmin, RHmax, RHmin, u2, Rn, G, P = x
    T = (Tmax + Tmin) / 2
    es = (e0(Tmax, xp) + e0(Tmin, xp)) / 2
    ea = (e0(Tmin, xp) * RHmax / 100 + e0(Tmax, xp) * RHmin / 100) / 2
    D = 4098 * e0(T, xp) / (T + 237.3) ** 2
    g = 0.665e-3 * P
    return (0.408 * D * (Rn - G) + g * 900 / (T + 273) * u2 * (es - ea)) / (D + g * (1 + 0.34 * u2))


def parede(K):                               # parede com K camadas: x = [Tin, Tout, hin, hout, L1..LK, k1..kK]
    def f(x, xp):
        R = 1 / x[2] + 1 / x[3]
        for Li, ki in zip(x[4:4 + K], x[4 + K:4 + 2 * K]):
            R = R + Li / ki
        return (x[0] - x[1]) / R
    return f


def monta_modelos(nomes):
    todos = {}
    for nome in nomes:
        if nome == "M1":
            todos[nome] = dict(f=M1, nomes=["P1", "V1", "V2", "cp", "cv"], r=[101.0, 0.1, 1.6, 1040.0, 740.0])
        elif nome == "M2":
            todos[nome] = dict(f=M2, nomes=["Tmax", "Tmin", "RHmax", "RHmin", "u2", "Rn", "G", "P"],
                               r=[30.0, 18.0, 85.0, 45.0, 2.0, 15.0, 0.6, 95.0])
        elif nome == "M3":
            todos[nome] = dict(f=parede(3), nomes=["Tin", "Tout", "hin", "hout", "L1", "L2", "L3", "k1", "k2", "k3"],
                               r=[40.0, 10.0, 8.0, 25.0, 0.02, 0.15, 0.05, 0.7, 1.0, 0.04])
        elif nome.startswith("W"):
            K = (int(nome[1:]) - 4) // 2
            todos[nome] = dict(f=parede(K), nomes=["Tin", "Tout", "hin", "hout"] + [f"L{i+1}" for i in range(K)] + [f"k{i+1}" for i in range(K)],
                               r=[40.0, 10.0, 8.0, 25.0] + [0.02 + 0.002 * i for i in range(K)] + [0.5 + 0.05 * i for i in range(K)])
        else:
            raise SystemExit(f"modelo desconhecido: {nome}")
    return todos


# ------------------------------------------------------------------ referência (60 dígitos) e medidas
def referencia(m):
    r = [mp.mpf(repr(v)) for v in m["r"]]
    return [mp.diff(lambda *a: m["f"](list(a), mp), r, tuple(1 if i == j else 0 for i in range(len(r)))) for j in range(len(r))]


def desvio(g, ref):
    return max(float(abs(mp.mpf(a) - b) / abs(b)) for a, b in zip(g, ref) if b != 0)


def melhor(fn, rep):
    ts, out = [], None
    for _ in range(rep):
        t = time.perf_counter(); out = fn(); ts.append(time.perf_counter() - t)
    return min(ts), out


# ------------------------------------------------------------------ métodos
def m_sympy(m, modo):
    X = sym.symbols(m["nomes"], real=True)
    fe = m["f"](list(X), sym)
    t0 = time.perf_counter(); grads = [sym.diff(fe, v) for v in X]; td = time.perf_counter() - t0
    val = dict(zip(X, m["r"]))
    if modo == "subs":
        te, g = melhor(lambda: [float(e.subs(val)) for e in grads], 3 if len(X) <= 24 else 1)
        return td, te, g
    t0 = time.perf_counter(); fl = sym.lambdify(X, grads, "math"); tl = time.perf_counter() - t0
    te, g = melhor(lambda: fl(*m["r"]), 200)
    return td + tl, te, list(g)


def m_symengine(m):
    X = sym.symbols(m["nomes"], real=True)
    sf = se.sympify(m["f"](list(X), sym)); SX = [se.Symbol(n) for n in m["nomes"]]
    t0 = time.perf_counter(); grads = [se.diff(sf, v) for v in SX]; td = time.perf_counter() - t0
    t0 = time.perf_counter(); lam = se.Lambdify(SX, grads); tl = time.perf_counter() - t0
    r = np.array(m["r"], dtype=float)
    te, g = melhor(lambda: lam(r), 200)
    return td + tl, te, [float(v) for v in np.ravel(g)]


def m_jax(m):
    gf = jax.jit(jax.grad(lambda x: m["f"](x, jnp)))
    r = jnp.array(m["r"], dtype=jnp.float64)
    t0 = time.perf_counter(); gf(r).block_until_ready(); tc = time.perf_counter() - t0
    te, g = melhor(lambda: gf(r).block_until_ready(), 300)
    return tc, te, [float(v) for v in g]


def m_autograd(m):
    gf = autograd.grad(lambda x: m["f"](x, anp))
    r = anp.array(m["r"], dtype=float)
    t0 = time.perf_counter(); gf(r); tc = time.perf_counter() - t0
    te, g = melhor(lambda: gf(r), 30)
    return tc, te, [float(v) for v in g]


def m_dif_finitas(m, modo="numpy"):
    """Diferenças centrais, h = 1e-5 * r.  modo="numpy": escalares NumPy (como no artigo);
    modo="float": ponto flutuante nativo do Python (3 a 7 vezes mais rápido)."""
    if modo == "numpy":
        r0 = np.array(m["r"], dtype=float); xp = np
        copia = lambda: r0.copy()
        lista = lambda v: list(v)
    else:
        r0 = list(m["r"]); xp = math
        copia = lambda: list(r0)
        lista = lambda v: v
    def grad():
        out = []
        for j in range(len(r0)):
            h = 1e-5 * abs(r0[j]); a = copia(); b = copia(); a[j] += h; b[j] -= h
            out.append((m["f"](lista(a), xp) - m["f"](lista(b), xp)) / (2 * h))
        return out
    te, g = melhor(grad, 5)
    return 0.0, te, g


def expr_maxima(nome, m):
    if nome == "M1":
        return "f: P1*(V1/V2)^(cp/cv)$"
    if nome == "M2":
        return ("e0(T):=0.6108*exp(17.27*T/(T+237.3))$ Tm:(Tmax+Tmin)/2$ es:(e0(Tmax)+e0(Tmin))/2$ "
                "ea:(e0(Tmin)*RHmax/100+e0(Tmax)*RHmin/100)/2$ Dl:4098*e0(Tm)/(Tm+237.3)^2$ gm:0.665e-3*P$ "
                "f:(0.408*Dl*(Rn-G)+gm*900/(Tm+273)*u2*(es-ea))/(Dl+gm*(1+0.34*u2))$")
    K = (len(m["nomes"]) - 4) // 2
    return "R: 1/hin+1/hout+" + "+".join(f"L{i+1}/k{i+1}" for i in range(K)) + "$ f:(Tin-Tout)/R$"


def m_maxima(nome, m):
    n = len(m["r"]); reps = 400 if n <= 14 else (100 if n <= 54 else 20)
    vals = ",".join(f"{a}={float(v)!r}" for a, v in zip(m["nomes"], m["r"]))
    roteiro = f"""fpprintprec:17$
{expr_maxima(nome, m)}
vl:[{",".join(m["nomes"])}]$
eqs:[{vals}]$
RR:{reps}$
t0:elapsed_real_time()$
for i:1 thru RR do grads:map(lambda([v],diff(f,v)),vl)$
t1:elapsed_real_time()$
for i:1 thru RR do nums:map(lambda([g],float(ev(g,eqs,numer))),grads)$
t2:elapsed_real_time()$
printf(true,"TD ~,6e~%",(t1-t0)/RR)$
printf(true,"TE ~,6e~%",(t2-t1)/RR)$
for v in nums do printf(true,"G ~,17e~%",v)$
"""
    with tempfile.NamedTemporaryFile("w", suffix=".mac", delete=False, encoding="utf-8") as fh:
        fh.write(roteiro); caminho = fh.name
    try:
        saida = subprocess.run([MAXIMA, "--very-quiet", "-b", caminho], capture_output=True, text=True, timeout=900).stdout
    finally:
        os.unlink(caminho)
    td = te = None; g = []
    for ln in saida.splitlines():
        if ln.startswith("TD "): td = float(ln[3:].replace("d", "e"))
        elif ln.startswith("TE "): te = float(ln[3:].replace("d", "e"))
        elif ln.startswith("G "): g.append(float(ln[2:].replace("d", "e")))
    if td is None or len(g) != n:
        raise RuntimeError("saída inesperada do Maxima")
    return td, te, g                                   # o Maxima já informa segundos por repetição


def metodos_disponiveis(max_subs):
    lista = [("SymPy (subs)", lambda nome, m: m_sympy(m, "subs") if len(m["r"]) <= max_subs else None),
             ("SymPy (lambdify)", lambda nome, m: m_sympy(m, "lam"))]
    if se is not None and np is not None: lista.append(("SymEngine", lambda nome, m: m_symengine(m)))
    if MAXIMA: lista.append(("Maxima", m_maxima))
    if jax is not None: lista.append(("JAX", lambda nome, m: m_jax(m)))
    if autograd is not None: lista.append(("Autograd", lambda nome, m: m_autograd(m)))
    if np is not None: lista.append(("Dif. finitas", lambda nome, m: m_dif_finitas(m, "numpy")))
    lista.append(("Dif. finitas (float)", lambda nome, m: m_dif_finitas(m, "float")))
    return lista


def ft(v):
    if v is None: return "—"
    s = f"{v:.0f}" if v >= 100 else f"{v:.1f}" if v >= 10 else f"{v:.2f}" if v >= 1 else f"{v:.3f}" if v >= 0.01 else f"{v:.4f}"
    return s.replace(".", ",")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--modelos", default="M1,M2,M3,W14,W24,W54")
    ap.add_argument("--max-subs", type=int, default=54)
    a = ap.parse_args()
    nomes = a.modelos.split(",")
    modelos = monta_modelos(nomes)
    metodos = metodos_disponiveis(a.max_subs)
    print("Métodos que serão executados:", ", ".join(n for n, _ in metodos))
    for rotulo, ativo in (("symengine", se), ("jax", jax), ("autograd", autograd), ("maxima", MAXIMA)):
        if not ativo: print(f"  (pulado: {rotulo} não encontrado)")
    # aquecimento: carrega módulos preguiçosos antes de medir
    aq = monta_modelos(["M1"])["M1"]
    for _, fn in metodos:
        try: fn("M1", aq)
        except Exception: pass
    res = {}
    for nome in nomes:
        m = modelos[nome]; ref = referencia(m); res[nome] = {"n": len(m["r"])}
        print(f"\nModelo {nome} (n = {len(m['r'])})")
        for rot, fn in metodos:
            try:
                o = fn(nome, m)
                if o is None:
                    print(f"  {rot:<21} não executado (n > --max-subs)"); continue
                prep, ev, g = o
                res[nome][rot] = dict(prep_ms=prep * 1e3, eval_ms=ev * 1e3, desvio=desvio(g, ref))
                print(f"  {rot:<21} prep {ft(prep*1e3):>8} ms | aval {ft(ev*1e3):>9} ms | desvio {res[nome][rot]['desvio']:.1e}")
            except Exception as e:
                print(f"  {rot:<21} falhou: {str(e)[:70]}")
    # tabelas finais
    rots = [n for n, _ in metodos]
    for titulo, chave, conv in (("Desvio relativo máximo", "desvio", lambda v: f"{v:.1e}"),
                                ("Tempo de preparação (ms)", "prep_ms", ft), ("Tempo de avaliação (ms)", "eval_ms", ft)):
        print(f"\n{titulo}\n{'Método':<22}" + "".join(f"{'n='+str(res[k]['n']):>11}" for k in nomes))
        for r_ in rots:
            print(f"{r_:<22}" + "".join(f"{(conv(res[k][r_][chave]) if r_ in res[k] else '—'):>11}" for k in nomes))
    res["_ambiente"] = dict(python=platform.python_version(), sistema=platform.platform(), sympy=sym.__version__,
                            mpmath=mp.__version__, symengine=getattr(se, "__version__", None), jax=getattr(jax, "__version__", None),
                            autograd=getattr(autograd, "__version__", None), maxima=bool(MAXIMA), nucleos=os.cpu_count())
    with open("resultados_desempenho.json", "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1, ensure_ascii=False)
    print("\nResultados gravados em resultados_desempenho.json")


if __name__ == "__main__":
    main()
