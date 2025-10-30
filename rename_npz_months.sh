#!/bin/bash

# Recursively rename NPZ files with three-digit month prefixes (001→01, 010→10, etc.)
# within the temporal directory tree.

ROOT="/datasets/work/oa-alantis/work/NESP_hydro/New_version/OCEANUS/temporal"

find "$ROOT" -type f -path "*/variables/*" -name '???_*_*_SS_Second_step.npz' -print0 |
while IFS= read -r -d '' f; do
  dir=$(dirname "$f")
  base=$(basename "$f")

  IFS=_ read -r m y var _rest <<< "$base"  # e.g., 001 2022 temp SS_Second_step.npz
  m2=$(printf "%02d" $((10#$m)))           # 001->01, 010->10, 012->12

  new="${m2}_${y}_${var}_SS_Second_step.npz"
  if [ "$base" != "$new" ]; then
    echo "Renaming: $base -> $new"
    mv -n -- "$f" "$dir/$new"
  fi
done


