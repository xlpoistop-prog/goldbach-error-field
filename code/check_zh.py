"""Cross-check that the Chinese version keeps every number of the English master."""
import re
import sys
import collections
import os

sys.stdout.reconfigure(encoding="utf-8")
PAPER_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "paper")

en = open(os.path.join(PAPER_DIR, "PAPER-EN.md"), encoding="utf-8").read()
zh = open(os.path.join(PAPER_DIR, "PAPER.md"), encoding="utf-8").read()

num = re.compile(r"[−-]?\d+(?:\.\d+)?(?:×10[⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+)?")

def toks(s):
    return collections.Counter(t.replace("−", "-") for t in num.findall(s))

a, b = toks(en), toks(zh)
missing = {k: v - b.get(k, 0) for k, v in a.items() if v > b.get(k, 0)}
extra = {k: v - a.get(k, 0) for k, v in b.items() if v > a.get(k, 0)}

# ignore trivial small integers that may legitimately be reworded (section refs, counts)
def interesting(t):
    return ("." in t) or ("×10" in t)

miss_i = {k: v for k, v in missing.items() if interesting(k)}
extra_i = {k: v for k, v in extra.items() if interesting(k)}

print("EN numeric tokens: %d distinct, %d total" % (len(a), sum(a.values())))
print("ZH numeric tokens: %d distinct, %d total" % (len(b), sum(b.values())))
print("\n-- numeric tokens in EN but missing/short in ZH (decimals only) --")
for k, v in sorted(miss_i.items(), key=lambda x: -x[1]):
    print("  %-14s missing %d of %d" % (k, v, a[k]))
print("\n-- numeric tokens in ZH not in EN (decimals only) --")
for k, v in sorted(extra_i.items(), key=lambda x: -x[1]):
    print("  %-14s extra %d" % (k, v))

# headings by numbering
def heads(s):
    out = []
    for ln in s.split("\n"):
        m = re.match(r"^(#{1,2})\s+(.*)$", ln)
        if m:
            num_ = re.match(r"^([0-9]+(?:\.[0-9]+)?|[A-C](?:\.[0-9]+)?)[\.\s]", m.group(2))
            out.append(num_.group(1) if num_ else m.group(2)[:30])
    return out

he, hz = heads(en), heads(zh)
print("\nEN headings: %d,  ZH headings: %d" % (len(he), len(hz)))
if he != hz:
    import difflib
    print("\n".join(list(difflib.unified_diff(he, hz, "EN", "ZH", lineterm=""))[:40]))
else:
    print("heading numbering identical")

# NOTE: the Chinese version is currently a digest, not a translation, so the
# decimal-token comparison is advisory. Flip this to an exit(1) once it is a
# full translation.
