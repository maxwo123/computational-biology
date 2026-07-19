"""
Google Cloud Batch task runner for Phase 6 docking.

Each Batch task receives BATCH_TASK_INDEX / BATCH_TASK_COUNT (injected by Batch)
and docks the stripe of the manifest assigned to it (round-robin by modulo, so
long and short ligands spread evenly across tasks). Inputs and outputs live in a
GCS bucket mounted at GCS_ROOT, so the pose/log cache works exactly like the
local pipeline: a task skips any (receptor, ligand) pair already computed.

Mirrors 04_scripts/06_dock_library.py — same vina invocation, same seed, same
pose/log naming (<receptor>__<ligand_id>.pdbqt / .log), so results are drop-in
compatible with 13_dock_candidates.py's cache and 14_rank_hits.py.
"""
import csv, os, subprocess, tempfile, time

ROOT = os.environ.get("GCS_ROOT", "/mnt/disks/gcs/hiv-docking")
IDX  = int(os.environ.get("BATCH_TASK_INDEX", "0"))
CNT  = int(os.environ.get("BATCH_TASK_COUNT", "1"))

# MANIFEST_FILE lets a triage job point at triage_manifest.csv instead of the
# default screen manifest (both live at the bucket root).
MANIFEST = os.path.join(ROOT, os.environ.get("MANIFEST_FILE", "manifest.csv"))
POSES    = os.path.join(ROOT, "poses")
os.makedirs(POSES, exist_ok=True)


def best_affinity(vina_stdout):
    for line in vina_stdout.splitlines():
        s = line.split()
        if len(s) >= 4 and s[0] == "1":
            return float(s[1])
    return None


def conf_for(receptor, exhaustiveness=None):
    """Rewrite the receptor path (and optionally exhaustiveness) into a temp conf."""
    src = os.path.join(ROOT, "confs", f"{receptor}_vina.conf")
    lines, saw_exh = [], False
    for l in open(src):
        low = l.strip().lower()
        if low.startswith("receptor"):
            lines.append(f"receptor = {os.path.join(ROOT, 'receptors', receptor + '.pdbqt')}\n")
        elif low.startswith("exhaustiveness") and exhaustiveness:
            lines.append(f"exhaustiveness = {exhaustiveness}\n"); saw_exh = True
        else:
            lines.append(l)
    if exhaustiveness and not saw_exh:
        lines.append(f"exhaustiveness = {exhaustiveness}\n")
    tf = tempfile.NamedTemporaryFile("w", suffix=".conf", delete=False)
    tf.writelines(lines); tf.close()
    return tf.name


rows = list(csv.DictReader(open(MANIFEST)))
mine = [r for i, r in enumerate(rows) if i % CNT == IDX]
print(f"[task {IDX}/{CNT}] {len(mine)} of {len(rows)} docking pairs")

t0 = time.time()
for n, r in enumerate(mine, 1):
    rec, lid = r["receptor"], r["ligand_id"]
    # optional per-row overrides (present in triage_manifest.csv, absent in the screen manifest)
    seed = (r.get("seed") or "1").strip()
    exh  = (r.get("exhaustiveness") or "").strip()
    tag  = f"__s{seed}" if r.get("seed") not in (None, "") else ""   # keep replicates distinct
    out = os.path.join(POSES, f"{rec}__{lid}{tag}.pdbqt")
    log = os.path.join(POSES, f"{rec}__{lid}{tag}.log")
    if os.path.exists(out) and os.path.exists(log):
        print(f"  [{n}/{len(mine)}] cached  {rec:13s} {lid}{tag}")
        continue
    lig = os.path.join(ROOT, "ligands", f"{lid}.pdbqt")
    if not os.path.exists(lig):
        print(f"  [{n}/{len(mine)}] MISSING ligand {lig} — skipping")
        continue
    conf = conf_for(rec, exh or None)
    v = subprocess.run(["vina", "--config", conf, "--ligand", lig,
                        "--out", out, "--seed", seed], capture_output=True, text=True)
    open(log, "w").write(v.stdout)
    aff = best_affinity(v.stdout)
    print(f"  [{n}/{len(mine)}] {rec:13s} {lid}{tag} {aff} kcal/mol")

print(f"[task {IDX}/{CNT}] done in {time.time()-t0:.0f}s")
