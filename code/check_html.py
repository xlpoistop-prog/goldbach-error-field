import os
import re
import sys
PAPER_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "paper")
sys.stdout.reconfigure(encoding="utf-8")
h = open(os.path.join(PAPER_DIR, "PAPER-EN.html"), encoding="utf-8").read()
print("imgs:", h.count("<img "), " figs:", h.count("<figure>"), " tables:", h.count("<table>"))
figs = [m.group(1) for m in re.finditer(r'<img src="([^"]+)"', h)]
for f in figs:
    print("  ", f)
bad = 0
for j, m in enumerate(re.finditer(r"<table>(.*?)</table>", h, re.S)):
    t = m.group(1)
    hdr = len(re.findall(r"<th>", t))
    rows = re.findall(r"<tr>(.*?)</tr>", t, re.S)
    counts = set(len(re.findall(r"<td>", r)) for r in rows if "<td>" in r)
    if counts != {hdr}:
        bad += 1
        print("BAD table", j, hdr, counts)
print("bad tables:", bad)

HTMLDIR = os.path.dirname(os.path.join(PAPER_DIR, "PAPER-EN.html"))
missing = [f for f in figs if not os.path.exists(os.path.join(HTMLDIR, f))]
if missing or bad:
    print("FAIL: missing figures %s, malformed tables %d" % (missing, bad))
    sys.exit(1)
print("OK: figures and tables check out")
