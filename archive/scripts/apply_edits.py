import ast, shutil, sys

SRC = "nutrient_tracer_fmt.py"
shutil.copy(SRC, SRC + ".bak")
lines = open(SRC).read().split("\n")

def find(sub, what):
    hits = [i for i, l in enumerate(lines) if sub in l]
    if len(hits) != 1:
        sys.exit("ABORT: %d candidates for %s, expected 1" % (len(hits), what))
    return hits[0]

# 1. colorbar tied to axes height
i = find("fig.colorbar(plt.cm.ScalarMappable", "colorbar call")
if "ax=ax" not in lines[i+1]:
    sys.exit("ABORT: unexpected line after colorbar: " + repr(lines[i+1]))
indent = lines[i][:len(lines[i]) - len(lines[i].lstrip())]
lines[i+1] = lines[i+1].replace("ax=ax", "cax=cax").replace(", pad=0.005", "")
lines.insert(i, indent + "cax = ax.inset_axes([1.02, 0, 0.015, 1], transform=ax.transAxes)")

# 2. two-line depth label (idempotent: handles edited or original text)
j = find("depth below surface (m)", "depth label")
lines[j] = lines[j].replace(
    ", profile stretched from AZMP 3-95 m",
    "\\n(AZMP 3-95 m compressed to 40 m)")
lines[j] = lines[j].replace(
    'depth below surface (m)"',
    'depth below surface (m)\\n(AZMP 3-95 m compressed to 40 m)"')

src = "\n".join(lines)
ast.parse(src)
open(SRC, "w").write(src)

print("OK. backup: " + SRC + ".bak")
print("--- colorbar ---")
for n in range(i-1, i+4):
    print("%d: %s" % (n+1, lines[n]))
print("--- label ---")
print("%d: %s" % (j+1, lines[j]))
