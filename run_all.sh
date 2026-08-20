#!/bin/sh
# Runs every verification script and writes one log per script into outputs/.
# Total run time is a few minutes, dominated by 05 (3132 parameter sets).
set -e
mkdir -p outputs figures
for f in scripts/*.py; do
    b=$(basename "$f" .py)
    printf '==> %s\n' "$b"
    python3 "$f" > "outputs/$b.log" 2>&1
done
cat outputs/*.log > outputs/ALL.log
printf 'done; logs in outputs/\n'
