def relatorio(fn, x, r, Δp, Δr, As, Rs, Pie, Me, Re, Mre, S, val):
    htm  = "<html><body><table border='1'>"

    htm += "<tr><td>x"
    for v in range(len(x)):
        htm += f"<td>{x[v]}"
    htm += "<td>Σ"

    htm += "<tr><td>r"
    for v in range(len(r)):
        htm += "<td>%.5f" % r[v]
    htm += "<td>"

    htm += "<tr><td>Δp"
    for v in range(len(Δp)):
        htm += "<td>%.5f" % Δp[v]
    htm += "<td>"

    htm += "<tr><td>Δr"
    for v in range(len(Δr)):
        htm += "<td>%.5f" % Δr[v]
    htm += "<td>"

    htm += "<tr><td>As"
    for v in range(len(As)):  htm += "<td>%.5f" % float(As[v].subs(val))
    htm += "<td>"

    htm += "<tr><td>Rs"
    for v in range(len(Rs)):  htm += "<td>%.5f" % float(Rs[v].subs(val))
    htm += "<td>%.5f" % float(S[0].subs(val))

    htm += "<tr><td>Pie"
    for v in range(len(Pie)): htm += "<td>%.5f" % float(Pie[v].subs(val))
    htm += "<td>%.5f" % float(S[1].subs(val))

    htm += "<tr><td>Me"
    for v in range(len(Me)):  htm += "<td>%.5f" % float(Me[v].subs(val))
    htm += "<td>%.5f" % float(S[2].subs(val))

    htm += "<tr><td>Re"
    for v in range(len(Re)):  htm += "<td>%.5f" % float(Re[v].subs(val))
    htm += "<td>%.5f" % float(S[3].subs(val))

    htm += "<tr><td>Mre"
    for v in range(len(Mre)): htm += "<td>%.5f" % float(Mre[v].subs(val))
    htm += "<td>%.5f" % float(S[4].subs(val))

    htm += "</table></body></html>"

    with open(fn, "w", encoding="utf-8") as f:
        f.write(htm.replace(".", ","))
