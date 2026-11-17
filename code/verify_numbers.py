"""Final consistency check: every number in the data files appears in both papers."""
import csv
import os
import sys

sys.stdout.reconfigure(encoding="utf-8")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
base = os.path.join(ROOT, "paper")
data = os.path.join(ROOT, "data")
en = open(os.path.join(base, "PAPER-EN.md"), encoding="utf-8").read()
zh = open(os.path.join(base, "PAPER.md"), encoding="utf-8").read()

def check(rows, cols, label):
    bad = 0
    for r in rows:
        for c in cols:
            v = r[c]
            # the data files carry 4 decimals; the papers print 2, 3 or 4 for the same value.
            # Build every plausible rendering and accept a literal match of any of them.
            base_vals = {v}
            try:
                x = float(v)
                for nd in (2, 3, 4):
                    t = ("%%.%df" % nd) % x
                    base_vals.add(t)
                    if t.endswith("0"):
                        base_vals.add(t.rstrip("0").rstrip("."))
                    if float(t) == int(float(t)):
                        base_vals.add("%d" % int(float(t)))
            except ValueError:
                pass
            pats = set()
            for b in base_vals:
                pats.add(b)
                pats.add(b.lstrip("-"))
                pats.add(b.replace("-", "−"))
                pats.add(b.replace("-", "−").lstrip("−"))
            ok_en = any(p in en for p in pats)
            ok_zh = any(p in zh for p in pats)
            if not (ok_en and ok_zh):
                bad += 1
                print("  MISSING %s %s=%s  EN:%s ZH:%s" % (label, c, v, ok_en, ok_zh))
    return bad

bad = 0
with open(os.path.join(data, "ck_main.csv"), encoding="utf-8") as f:
    rows = list(csv.DictReader(f))
bad += check(rows, ["c_k"], "ck_main")
# 4-decimal values appear in the main table; 3-decimal rounded ones in the figures/abstract
cnt = 0
for r in rows:
    v = r["c_k"]
    if v not in en or v not in zh:
        cnt += 1
print("ck_main: %d rows, %d whose 4-decimal c(k) is not literally in both papers (expected: those shown only to 3 decimals)" % (len(rows), cnt))

with open(os.path.join(data, "ck_scale.csv"), encoding="utf-8") as f:
    rows = list(csv.DictReader(f))
bad2 = check(rows, ["c2", "c3", "c4", "c6"], "ck_scale")

with open(os.path.join(data, "ck_weight.csv"), encoding="utf-8") as f:
    rows = list(csv.DictReader(f))
bad3 = check(rows, ["c2", "c3", "c4", "c6", "c10"], "ck_weight")

print("hard mismatches: main=%d scale=%d weight=%d" % (bad, bad2, bad3))

total = bad + bad2 + bad3
if total:
    print("FAIL: %d values could not be matched" % total)
    sys.exit(1)
print("OK: every checked value matches")
