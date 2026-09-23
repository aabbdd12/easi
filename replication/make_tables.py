# -*- coding: utf-8 -*-
r"""Build the LaTeX tables of the technical note from the e() matrices that
tests/../docs/fig/*.csv hold (written by the Mata helper in the note's
do-file).  Run:  python make_tables.py
"""
import csv
import io
import os
import sys

# directory holding the elast_*.csv written by 11_hixdata_results.do; the
# tab_*.tex files are written next to them.  Default: this script's folder.
d = (sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))) + os.sep


def rd(f):
    rows = list(csv.reader(io.open(d + f, encoding="utf-8")))
    cols = rows[0][1:]
    rn = [r[0] for r in rows[1:]]
    M = [[float(x) for x in r[1:]] for r in rows[1:]]
    return rn, cols, M


def esc(s):
    return s.replace("_", r"\_").replace("#", r"\#")


def mat_table(name, fe, fs, caption, label, rowhdr, colhdr, fmt="%7.3f"):
    rn, cn, E = rd(fe)
    _, _, S = rd(fs)
    n = len(cn)
    out = [r"\begin{table}[ht]", r"\centering", r"\footnotesize",
           r"\caption{%s}" % caption, r"\label{%s}" % label,
           r"\setlength{\tabcolsep}{3.5pt}",
           r"\begin{tabular}{l" + "r" * n + "}", r"\toprule",
           rowhdr + r" & \multicolumn{%d}{c}{%s} \\" % (n, colhdr),
           r"\cmidrule(lr){2-%d}" % (n + 1),
           " & " + " & ".join(esc(c) for c in cn) + r" \\", r"\midrule"]
    for i, r in enumerate(rn):
        out.append(esc(r) + " & " + " & ".join(fmt % v for v in E[i]) + r" \\")
        out.append(" & " + " & ".join(r"{\scriptsize(%s)}" % (fmt % v).strip() for v in S[i]) + r" \\[1pt]")
    out += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    io.open(d + name, "w", encoding="utf-8").write("\n".join(out))
    print(name, len(rn), "x", n)


mat_table("tab_price_nc.tex", "elast_price_nc.csv", "elast_price_nc_se.csv",
          r"Uncompensated (Marshallian) price elasticities, $\partial\ln q_j/\partial\ln p_k$, "
          r"with delta-method standard errors in parentheses. \emph{hixdata}, reference specification.",
          "tab:out_price_nc", "good", "price of")
mat_table("tab_price_c.tex", "elast_price_c.csv", "elast_price_c_se.csv",
          r"Compensated (Hicksian) price elasticities at constant utility, with standard errors. "
          r"\emph{hixdata}, reference specification.",
          "tab:out_price_c", "good", "price of")
mat_table("tab_demo.tex", "elast_demo.csv", "elast_demo_se.csv",
          r"Demographic elasticities $\partial w_j/\partial z_t$ (semi-elasticities of the shares), "
          r"with standard errors. \emph{hixdata}, reference specification.",
          "tab:out_demo", "demographic", "good", fmt="%8.4f")

rn, cn, E = rd("elast_exp.csv")
_, _, S = rd("elast_exp_se.csv")
out = [r"\begin{table}[ht]", r"\centering",
       r"\caption{Expenditure elasticities $\partial\ln q_j/\partial\ln x$ with delta-method standard "
       r"errors. \emph{hixdata}, reference specification; the last good is recovered by adding up.}",
       r"\label{tab:out_exp}", r"\begin{tabular}{lrr}", r"\toprule",
       r"good & elasticity & std.\ err. \\", r"\midrule"]
for j, c in enumerate(cn):
    out.append("%s & %.4f & %.4f " % (esc(c), E[0][j], S[0][j]) + r"\\")
out += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
io.open(d + "tab_exp.tex", "w", encoding="utf-8").write("\n".join(out))
print("tab_exp.tex", len(cn))
