# Phase 6 docking on Google Cloud Batch (Spot VMs)

Fans out the ~285 (95 candidates × 3 receptors) AutoDock Vina dockings across
many cheap **Spot** VMs, finishing in **~15–25 min** instead of most of a day —
using the **same Vina 1.2.7** as the local pipeline, so scores stay identical.

Only the **docking** step runs in the cloud. You still build/filter/prep the
library locally (steps 10–12, they need RDKit/Meeko + PubChem) and rank locally
(step 14). The cloud only replaces step 13.

## Why this design
- **Embarrassingly parallel:** every (receptor, ligand) docking is independent.
- **Spot-safe:** each pose/log is cached in the bucket, so if a Spot VM is
  preempted the job just re-runs the unfinished pairs — no lost work.
- **No new science risk:** identical Vina binary + seed → comparable to your
  local ritonavir 1.93 Å benchmark and Phase 4 numbers. (A GPU engine like
  Uni-Dock would be faster for *huge* libraries but changes scores → re-validate.)

## Prerequisites
- A GCP project with **billing enabled**.
- Run everything from **Google Cloud Shell** (https://shell.cloud.google.com) —
  it has `gcloud`, `gsutil`, `envsubst`, and Docker preinstalled, so you don't
  install anything locally. Upload this repo (or just the `gcp/`, `01_proteins/`,
  `02_ligands/candidates/`, `03_docking/*.conf` files) to Cloud Shell.
- To pull results back to your Mac afterwards, install the SDK locally once:
  `brew install --cask google-cloud-sdk`  (gives you `gcloud storage`).

## Files
| file | role |
|---|---|
| `Dockerfile` | tiny image: Vina 1.2.7 + the task runner |
| `batch_dock_task.py` | one Batch task = docks its stripe of the manifest, caches to the bucket |
| `make_manifest.py` | builds `manifest.csv` (receptor × prepared-candidate pairs) |
| `batch_job.json.template` | Batch job spec (Spot, GCS mount) — filled by `envsubst` |
| `run_on_gcp.sh` | orchestrator: `setup → build → stage → submit → status → fetch` |

## Run it
```bash
cd gcp
export PROJECT_ID=your-project BUCKET=your-unique-bucket REGION=us-central1
# optional tuning:  export TASK_COUNT=50 PARALLELISM=50   # PARALLELISM×4 = concurrent vCPUs

./run_on_gcp.sh setup     # enable APIs, create Artifact Registry repo + bucket
./run_on_gcp.sh build     # Cloud Build the image (a few min)
./run_on_gcp.sh stage     # upload receptors, ligands, confs, manifest
./run_on_gcp.sh submit    # launch the Batch job
./run_on_gcp.sh status    # repeat until state = SUCCEEDED
./run_on_gcp.sh fetch     # download poses into ../03_docking/poses/
```
Then locally, turn the poses back into results:
```bash
.venv/bin/python 04_scripts/13_dock_candidates.py   # all cached → rebuilds docking_results_candidates.csv instantly
.venv/bin/python 04_scripts/14_rank_hits.py         # ranking + fig5
```

## Cost & quota notes
- `c2d-standard-4` Spot is typically well under **$0.05/hr**; 50 of them for ~20 min
  is roughly **a dollar or two** for the whole screen. Verify current prices in the
  console — Spot rates fluctuate.
- `PARALLELISM=50` × 4 vCPU = **200 concurrent vCPUs**. New projects often have a
  lower regional Spot vCPU quota — if `submit` complains, lower `PARALLELISM` or
  request a quota bump (IAM & Admin → Quotas).
- Delete the bucket when done to avoid storage charges: `gcloud storage rm -r gs://$BUCKET`.

## Top-hit triage on the same setup (Step 32)
After the screen, `04_scripts/15_triage_hits.py` re-docks the best hits with
several random seeds at high exhaustiveness to keep only *robust* predictions.
It emits `gcp/triage_manifest.csv` (columns `receptor, ligand_id, seed,
exhaustiveness`), which the same `batch_dock_task.py` runs unchanged — just point
the job at it with `MANIFEST_FILE`:

```bash
# 1) locally, after the screen's step 14 has written candidate_hits.csv:
.venv/bin/python 04_scripts/15_triage_hits.py manifest      # -> gcp/triage_manifest.csv

# 2) on Batch (ligands already staged from the screen; this just adds the manifest):
cd gcp
MANIFEST_FILE=triage_manifest.csv JOB_NAME=hiv-phase6-triage TASK_COUNT=30 PARALLELISM=30 ./run_on_gcp.sh stage
MANIFEST_FILE=triage_manifest.csv JOB_NAME=hiv-phase6-triage TASK_COUNT=30 PARALLELISM=30 ./run_on_gcp.sh submit
./run_on_gcp.sh fetch

# 3) locally, turn replicate poses into confidence scores + fig6:
.venv/bin/python 04_scripts/15_triage_hits.py aggregate     # -> 05_results/triage_confidence.csv
```
Seed-tagged poses (`<rec>__<lid>__s<seed>.pdbqt`) don't collide with the screen's
poses, so the cache stays clean. No cloud? `15_triage_hits.py dock-local` runs the
identical dockings on your Mac.

## Scaling further
This exact setup handles 300 or 300,000 dockings — bump `taskCount`/`PARALLELISM`.
