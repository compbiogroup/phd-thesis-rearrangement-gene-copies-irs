#!/usr/bin/env bash
# Convert every page of this folder to PDF (in pdf/), using pandoc and weasyprint.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p pdf
for f in *.md; do
  n=${f%.md}
  pandoc "$f" -s -o "pdf/$n.pdf" --pdf-engine=weasyprint \
    --lua-filter=pdf-build/links.lua --css=pdf-build/doc.css --metadata pagetitle="$n"
done
