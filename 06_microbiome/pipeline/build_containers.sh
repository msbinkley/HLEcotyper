#!/bin/bash
# build_containers.sh -- build pinned Kraken2 2.17.1 + Bracken 3.1 Apptainer image.
#
# WHERE TO RUN: a node WITH internet access. On Sherlock that means a LOGIN node or
# the sh-dtn data-transfer node -- compute nodes have NO internet, so `apptainer build`
# (which pulls the ubuntu base + downloads release tarballs) will fail there.
# Use tmux/screen so it survives disconnects.
set -euo pipefail

TOOLS=/oak/stanford/groups/emoding/analysis/brian/tools/microbiome
SIF="$TOOLS/containers/kraken2_bracken_2.17.1_3.1.sif"
DEF="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/kraken2_bracken.def"

mkdir -p "$TOOLS/containers"
# Keep Apptainer's cache + scratch OFF your small $HOME quota.
export APPTAINER_CACHEDIR="$TOOLS/.apptainer/cache"
export APPTAINER_TMPDIR="$TOOLS/.apptainer/tmp"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"

# Load the container runtime (adjust to your cluster's module name).
module load system apptainer 2>/dev/null || module load apptainer 2>/dev/null || ml apptainer || true

echo ">>> Building $SIF from $DEF"
# --fakeroot lets an unprivileged user run the %post apt-get/compile steps.
apptainer build --fakeroot "$SIF" "$DEF"

echo ">>> Verifying versions inside the image:"
apptainer exec "$SIF" kraken2 --version | head -1
apptainer exec "$SIF" bracken -v

cat <<EOF

DONE. New combined image:
  $SIF

Point BOTH variables in config.sh at it (one image provides both tools):
  export KRAKEN2_SIF="$SIF"
  export BRACKEN_SIF="$SIF"

------------------------------------------------------------------------------
PLAN B -- if your cluster does NOT permit '--fakeroot' / unprivileged .def builds,
pull prebuilt images instead (tags may lag the exact upstream version; check the
registry's tag list and pick the newest 2.17.x / 3.1):
  apptainer build "$TOOLS/containers/kraken2.sif" docker://staphb/kraken2:latest
  apptainer build "$TOOLS/containers/bracken.sif" docker://staphb/bracken:latest
  # then verify: apptainer exec .../kraken2.sif kraken2 --version
  #              apptainer exec .../bracken.sif bracken -v
  # BioContainers alternative: docker://quay.io/biocontainers/kraken2:<tag>
  #                            docker://quay.io/biocontainers/bracken:<tag>
------------------------------------------------------------------------------
EOF
