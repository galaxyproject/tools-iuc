#!/usr/bin/env bash
set -euo pipefail
# Assume we have the original substrate_dir.zip with numbered files
unzip -q substrate_dir.zip -d substrate_tmp
cd substrate_tmp
cp glucose_1.sdf glucose.sdf
cp fructose_1.sdf fructose.sdf
zip -q -r ../substrate_dir.zip glucose.sdf fructose.sdf
cd ..
rm -rf substrate_tmp
