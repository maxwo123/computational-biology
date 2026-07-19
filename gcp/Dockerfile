# Minimal image for GPU-free AutoDock Vina docking tasks on Google Cloud Batch.
# Contains ONLY the Vina binary + the task runner — receptors/ligands live in a
# mounted GCS bucket, not baked into the image.
FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        wget ca-certificates libstdc++6 \
    && rm -rf /var/lib/apt/lists/*

# AutoDock Vina 1.2.7 static Linux binary (same version validated locally).
ARG VINA_VERSION=1.2.7
RUN wget -q "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v${VINA_VERSION}/vina_${VINA_VERSION}_linux_x86_64" \
        -O /usr/local/bin/vina \
    && chmod +x /usr/local/bin/vina \
    && vina --version

COPY batch_dock_task.py /app/batch_dock_task.py
WORKDIR /app
ENTRYPOINT ["python3", "/app/batch_dock_task.py"]
