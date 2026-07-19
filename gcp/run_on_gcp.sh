#!/usr/bin/env bash
# Orchestrates the Phase 6 docking fan-out on Google Cloud Batch (Spot VMs).
#
# Run the subcommands in order:  setup -> build -> stage -> submit -> status -> fetch
# All GCP-specific values come from the env vars below — edit them or export first.
#
#   cd gcp
#   export PROJECT_ID=my-project BUCKET=my-bucket REGION=us-central1
#   ./run_on_gcp.sh setup     # enable APIs + create Artifact Registry repo + bucket
#   ./run_on_gcp.sh build     # Cloud Build the Vina image (no local Docker needed)
#   ./run_on_gcp.sh stage     # upload receptors/ligands/confs/manifest to the bucket
#   ./run_on_gcp.sh submit    # launch the Batch job
#   ./run_on_gcp.sh status    # watch it
#   ./run_on_gcp.sh fetch     # pull poses back into ../03_docking/poses/
set -euo pipefail

PROJECT_ID="${PROJECT_ID:?set PROJECT_ID}"
BUCKET="${BUCKET:?set BUCKET (name only, no gs://)}"
REGION="${REGION:-us-central1}"
REPO="${REPO:-hiv-docking}"
JOB_NAME="${JOB_NAME:-hiv-phase6-dock}"
TASK_COUNT="${TASK_COUNT:-50}"      # number of Batch tasks (manifest is striped across them)
PARALLELISM="${PARALLELISM:-50}"   # how many run at once (watch your Spot vCPU quota: x4 vCPU each)
MANIFEST_FILE="${MANIFEST_FILE:-manifest.csv}"  # set to triage_manifest.csv for the triage job

HERE="$(cd "$(dirname "$0")" && pwd)"
PROJ="$(cd "$HERE/.." && pwd)"
GPREFIX="gs://${BUCKET}/hiv-docking"

case "${1:-}" in
setup)
  gcloud config set project "$PROJECT_ID"
  gcloud services enable batch.googleapis.com compute.googleapis.com \
      artifactregistry.googleapis.com cloudbuild.googleapis.com storage.googleapis.com
  gcloud artifacts repositories create "$REPO" --repository-format=docker \
      --location="$REGION" 2>/dev/null || echo "repo exists"
  gcloud storage buckets create "gs://${BUCKET}" --location="$REGION" 2>/dev/null || echo "bucket exists"
  ;;
build)
  gcloud builds submit "$HERE" \
      --tag "${REGION}-docker.pkg.dev/${PROJECT_ID}/${REPO}/hiv-vina:latest"
  ;;
stage)
  # (re)build the manifest from whatever candidates are prepared locally.
  # make_manifest.py is stdlib-only, so fall back to system python3 in Cloud
  # Shell, where the project venv does not exist.
  PY="${PROJ}/.venv/bin/python"
  [ -x "$PY" ] || PY="$(command -v python3)"
  "$PY" "${HERE}/make_manifest.py"
  gcloud storage cp "${PROJ}/01_proteins/hiv_protease.pdbqt" \
                    "${PROJ}/01_proteins/cathepsin_d.pdbqt" \
                    "${PROJ}/01_proteins/bace1.pdbqt"          "${GPREFIX}/receptors/"
  gcloud storage cp "${PROJ}/03_docking/"*_vina.conf           "${GPREFIX}/confs/"
  gcloud storage cp "${PROJ}/02_ligands/candidates/"*.pdbqt    "${GPREFIX}/ligands/"
  gcloud storage cp "${HERE}/manifest.csv"                     "${GPREFIX}/manifest.csv"
  [ -f "${HERE}/triage_manifest.csv" ] && \
    gcloud storage cp "${HERE}/triage_manifest.csv"            "${GPREFIX}/triage_manifest.csv" || true
  echo "staged inputs to ${GPREFIX}/"
  ;;
submit)
  export PROJECT_ID BUCKET REGION REPO TASK_COUNT PARALLELISM MANIFEST_FILE
  envsubst < "${HERE}/batch_job.json.template" > "${HERE}/batch_job.json"
  gcloud batch jobs submit "$JOB_NAME" --location "$REGION" \
      --config "${HERE}/batch_job.json"
  echo "submitted ${JOB_NAME}"
  ;;
status)
  gcloud batch jobs describe "$JOB_NAME" --location "$REGION" \
      --format="value(status.state, status.taskGroups)"
  ;;
fetch)
  gcloud storage cp "${GPREFIX}/poses/*" "${PROJ}/03_docking/poses/"
  echo "poses pulled -> ${PROJ}/03_docking/poses/"
  echo "now locally:  .venv/bin/python 04_scripts/13_dock_candidates.py  (rebuilds CSV from cache)"
  echo "then:         .venv/bin/python 04_scripts/14_rank_hits.py"
  ;;
*)
  echo "usage: $0 {setup|build|stage|submit|status|fetch}"; exit 1;;
esac
