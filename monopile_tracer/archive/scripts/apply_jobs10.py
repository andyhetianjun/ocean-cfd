import ast, shutil, sys
SRC = "nutrient_tracer_fmt.py"
src = open(SRC).read()

old_lab = '    lab = UNITS[NUTRIENT]; nut = NUTRIENT.capitalize()'
new_lab = ('    lab = "Nutrients (mmol m$^{-3}$)"; nut = "Nutrients"; fp = "nutrients"\n'
           '    # NUTRIENT stays "nitrate" for the AZMP column lookup; only display/filenames change')
assert src.count(old_lab) == 1, "lab line not unique"
src = src.replace(old_lab, new_lab)

old_jobs = '''    jobs = [("xz", 60.0, f"{NUTRIENT}_xz_y60_centre.gif",
             f"{nut}  |  vertical X-Z slice at y = 60 m (centre, through monopile)"),
            ("xz", 52.0, f"{NUTRIENT}_xz_y52_minusD.gif",
             f"{nut}  |  vertical X-Z slice at y = 52 m (centre - D)"),
            ("xz", 68.0, f"{NUTRIENT}_xz_y68_plusD.gif",
             f"{nut}  |  vertical X-Z slice at y = 68 m (centre + D)"),
            ("xy", float(z_tr[-1]), f"{NUTRIENT}_xy_top.gif",
             f"{nut}  |  horizontal X-Y plane at top Z (z = +20 m, surface)"),
            ("xy", float(z_tr[len(z_tr)//2]), f"{NUTRIENT}_xy_middle.gif",
             f"{nut}  |  horizontal X-Y plane at mid Z (z = 0 m)"),
            ("xy", float(z_tr[0]), f"{NUTRIENT}_xy_bottom.gif",
             f"{nut}  |  horizontal X-Y plane at bottom Z (z = -20 m)")]'''

new_jobs = '''    jobs = [("xz", 56.0, f"{fp}_xz_y56.gif",
             f"{nut}  |  vertical X-Z slice at y = 56 m (centre - 4 m, pile edge)"),
            ("xz", 58.0, f"{fp}_xz_y58.gif",
             f"{nut}  |  vertical X-Z slice at y = 58 m (centre - 2 m)"),
            ("xz", 60.0, f"{fp}_xz_y60.gif",
             f"{nut}  |  vertical X-Z slice at y = 60 m (centre, through monopile)"),
            ("xz", 62.0, f"{fp}_xz_y62.gif",
             f"{nut}  |  vertical X-Z slice at y = 62 m (centre + 2 m)"),
            ("xz", 64.0, f"{fp}_xz_y64.gif",
             f"{nut}  |  vertical X-Z slice at y = 64 m (centre + 4 m, pile edge)"),
            ("xy", float(z_tr[-1]), f"{fp}_xy_d00_surface.gif",
             f"{nut}  |  horizontal X-Y plane at 0 m depth (surface)"),
            ("xy", 10.0, f"{fp}_xy_d10_subsurface.gif",
             f"{nut}  |  horizontal X-Y plane at 10 m depth (subsurface)"),
            ("xy", 0.0, f"{fp}_xy_d20_middle.gif",
             f"{nut}  |  horizontal X-Y plane at 20 m depth (middle)"),
            ("xy", -10.0, f"{fp}_xy_d30_subbottom.gif",
             f"{nut}  |  horizontal X-Y plane at 30 m depth (sub-bottom)"),
            ("xy", float(z_tr[0]), f"{fp}_xy_d40_bottom.gif",
             f"{nut}  |  horizontal X-Y plane at 40 m depth (bottom)")]'''

assert src.count(old_jobs) == 1, "jobs block not matched verbatim"
src = src.replace(old_jobs, new_jobs)

ast.parse(src)
shutil.copy(SRC, SRC + ".bak_jobs")
open(SRC, "w").write(src)
print("OK, backup nutrient_tracer_fmt.py.bak_jobs")
