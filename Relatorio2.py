def relatorio2(fn, x, r, Δp, Δr, saidas):
    f   = lambda v: "%.5f" % v
    tds = lambda v: "".join("<td>" + f(e) for e in v)
    linha = lambda nome, v, s="": "<tr><td>" + nome + tds(v) + "<td>" + s
    htm = "<html><body>"

    for s in saidas:
        htm += "<h3>" + s["nome"] + " = " + f(s["f0"]) + "</h3><table border='1'>"
        htm += "<tr><td>x" + "".join(f"<td>{v}" for v in x) + "<td>Σ"
        htm += linha("r",  r)
        htm += linha("Δp", Δp)
        htm += linha("Δr", Δr)
        htm += linha("As", s["As"])
        htm += linha("Rs", s["Rs"], f(sum(s["Rs"])))
        htm += linha("Pie", s["Pie"], f(sum(s["Pie"])))
        htm += linha("Me",  s["Me"],  f(sum(s["Me"])))
        htm += linha("Re",  s["Re"],  f(sum(s["Re"])))
        htm += linha("Mre", s["Mre"], f(sum(s["Mre"])))
        htm += "</table>"

        p1, m1 = sum(s["Pie"]), sum(s["Me"])
        p2, m2, f0 = s["Pie2"], s["Me2"], s["f0"]
        htm += "<table border='1'><tr><td>ordem<td>Pie<td>Me<td>Re<td>Mre"
        htm += "<tr><td>1" + tds((p1, m1, p1/f0, m1/f0))
        htm += "<tr><td>2" + tds((p2, m2, p2/f0, m2/f0))
        htm += "</table>"

    htm += "</body></html>"

    with open(fn, "w", encoding="utf-8") as g:
        g.write(htm.replace(".", ","))
