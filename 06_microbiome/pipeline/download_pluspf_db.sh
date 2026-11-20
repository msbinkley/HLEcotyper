#!/bin/bash
# download_pluspf_db.sh -- fetch + verify + extract the latest PlusPF Kraken2/Bracken DB
# (Standard + fungi + protozoa), which the bundled Bracken k-mer distributions make
# ready for READ_LEN=100 with no rebuild.
#
# WHERE TO RUN: a node WITH internet (Sherlock LOGIN node or sh-dtn). ~80 GB download +
# ~150 GB extracted -- use tmux/screen. Compute nodes have no internet.
set -euo pipefail

# >>> Set DATE to the NEWEST build listed at https://benlangmead.github.io/aws-indexes/k2
#     (builds are quarterly; 20260226 was latest as of Feb 2026 -- check for newer).
DATE="20260226"

TOOLS=/oak/stanford/groups/emoding/analysis/brian/tools/microbiome
DEST="$TOOLS/db/k2_pluspf_${DATE}"
TARBALL="k2_pluspf_${DATE}.tar.gz"
URL="https://genome-idx.s3.amazonaws.com/kraken/${TARBALL}"

mkdir -p "$DEST"
cd "$DEST"

echo ">>> Downloading $URL  (resumable; ~80 GB)"
wget -c "$URL"

echo ">>> Checksum (compare this to the MD5 shown on the genome-idx web page):"
md5sum "$TARBALL"

echo ">>> Extracting (files unpack loose into $DEST)"
tar -xzvf "$TARBALL"

echo ">>> Sanity check -- required index + Bracken 100-mer distribution present:"
ls -lh hash.k2d taxo.k2d opts.k2d database100mers.kmer_distrib inspect.txt

echo ">>> Library composition (should now include fungi + protozoa):"
if [ -f library_report.tsv ]; then cut -f1 library_report.tsv | sort | uniq -c; fi

cat <<EOF

DONE. New database directory:
  $DEST

Update config.sh:
  export KRAKEN2_DB="$DEST"

Optional: once verified, delete the tarball to reclaim ~80 GB:
  rm -f "$DEST/$TARBALL"
EOF
