import ast, shutil, sys
SRC = "nutrient_tracer_fmt.py"
src = open(SRC).read()

def sub(old, new, what):
    global src
    n = src.count(old)
    if n != 1: sys.exit("ABORT: %s found %d times, expected 1" % (what, n))
    src = src.replace(old, new)

sub('tt.set_text(f"t = {int(i*DT_S)} s")',
    'tt.set_text(f"t = {int(TIMESTEPS[0]) + int(i*DT_S)} s")', 'timestamp')

sub('''    import time
    t0 = time.time()
    print(f"advecting {len(TIMESTEPS)} steps on fine tracer grid...", flush=True)
    all_slabs = evolve(C0, x, y, z_vel, z_tr, [(m, sv, fn) for m, sv, fn, _ in jobs])
    print(f"  done in {time.time()-t0:.0f}s", flush=True)''',
'''    import time
    CACHE = "slabs_cache"
    sig = "|".join([VEL_DIR, str(TRACER_NZ), TIMESTEPS[0], TIMESTEPS[-1],
                    str(len(TIMESTEPS))] + [f"{m}:{sv}" for m, sv, _, _ in jobs])
    sigf = os.path.join(CACHE, "sig.txt")
    all_slabs = None
    if os.path.exists(sigf) and open(sigf).read() == sig:
        all_slabs = [list(np.load(os.path.join(CACHE, f"arr_{k}.npy")))
                     for k in range(len(jobs))]
        print(f"  loaded cached slabs from {CACHE}/ (advection skipped)", flush=True)
    elif os.path.exists(sigf):
        print(f"  {CACHE}/ exists but config changed, re-advecting", flush=True)
    if all_slabs is None:
        t0 = time.time()
        print(f"advecting {len(TIMESTEPS)} steps on fine tracer grid...", flush=True)
        all_slabs = evolve(C0, x, y, z_vel, z_tr, [(m, sv, fn) for m, sv, fn, _ in jobs])
        print(f"  done in {time.time()-t0:.0f}s", flush=True)
        t0 = time.time()
        os.makedirs(CACHE, exist_ok=True)
        for k, s in enumerate(all_slabs):
            np.save(os.path.join(CACHE, f"arr_{k}.npy"), np.stack(s))
        open(sigf, "w").write(sig)
        print(f"  cached slabs to {CACHE}/ ({time.time()-t0:.0f}s)", flush=True)''', 'slab cache')

ast.parse(src)
shutil.copy(SRC, SRC + ".bak_cache")
open(SRC, "w").write(src)
print("OK, backup nutrient_tracer_fmt.py.bak_cache")
