# -*- coding: utf-8 -*-
"""README.txt from README.md: pipe tables become indented lists (pandoc's
plain writer keeps them as 250-character rows), everything else goes
through pandoc.  Run from the repository root:  python replication/make_readme_txt.py"""
import io, re, subprocess, sys, os

def table_to_list(rows):
    cells = [[c.strip() for c in r.strip().strip("|").split("|")] for r in rows]
    header, body = cells[0], [c for c in cells[1:] if not set("".join(c)) <= set("-: ")]
    out = []
    for r in body:
        out.append("- " + r[0])
        for h, c in zip(header[1:], r[1:]):
            if c and c != "—":
                out.append("    - " + h + ": " + c)
        out.append("")
    return out

def convert(md):
    lines = md.splitlines(); out = []; i = 0
    while i < len(lines):
        if lines[i].lstrip().startswith("|"):
            j = i
            while j < len(lines) and lines[j].lstrip().startswith("|"): j += 1
            out += table_to_list(lines[i:j]); i = j
        else:
            out.append(lines[i]); i += 1
    return "\n".join(out)

for d in [".", "replication"]:
    md = io.open(os.path.join(d, "README.md"), encoding="utf-8").read()
    txt = subprocess.run(["pandoc", "-f", "gfm", "-t", "plain", "--columns=78"],
                         input=convert(md).encode("utf-8"), capture_output=True, check=True).stdout.decode("utf-8")
    txt = "(generated from README.md; edit that file)\n\n" + txt
    io.open(os.path.join(d, "README.txt"), "w", encoding="utf-8", newline="\n").write(txt)
    print(d + "/README.txt", len(txt.splitlines()), "lines")
