#!/usr/bin/env bash
set -euo pipefail
base="https://raw.githubusercontent.com/bioinfbrad/enzywizard-mut-integrate/main/examples"
rm -rf mut_integrate_wt_input mut_integrate_mut_input
mkdir -p mut_integrate_wt_input mut_integrate_mut_input
wt_files=(
  aaprops_report_cleaned_1ZG4_WT.json
  conservation_report_cleaned_1ZG4_WT.json
  disorder_report_cleaned_1ZG4_WT.json
  dock_report_cleaned_1ZG4_WT_glucose_fructose.json
  embedding_report_cleaned_1ZG4_WT.json
  energy_report_cleaned_1ZG4_WT.json
  flexibility_report_cleaned_1ZG4_WT.json
  hydrocluster_report_cleaned_1ZG4_WT.json
  interaction_report_cleaned_1ZG4_WT_docked_glucose_docked_fructose.json
  pocket_report_cleaned_1ZG4_WT.json
  substrate_report_glucose_fructose.json
)
mut_files=(
  aaprops_report_cleaned_1ZG6_S70G.json
  conservation_report_cleaned_1ZG6_S70G.json
  disorder_report_cleaned_1ZG6_S70G.json
  dock_report_cleaned_1ZG6_S70G_glucose_fructose.json
  embedding_report_cleaned_1ZG6_S70G.json
  energy_report_cleaned_1ZG6_S70G.json
  flexibility_report_cleaned_1ZG6_S70G.json
  hydrocluster_report_cleaned_1ZG6_S70G.json
  interaction_report_cleaned_1ZG6_S70G_docked_glucose_docked_fructose.json
  pocket_report_cleaned_1ZG6_S70G.json
  substrate_report_glucose_fructose.json
)
for f in "${wt_files[@]}"; do curl -fsSL "$base/wt_input/$f" -o "mut_integrate_wt_input/$f"; done
for f in "${mut_files[@]}"; do curl -fsSL "$base/mut_input/$f" -o "mut_integrate_mut_input/$f"; done
rm -f wt_input.zip mut_input.zip
find mut_integrate_wt_input mut_integrate_mut_input -type f -name 'mut_clean_report*.json' -delete
(cd mut_integrate_wt_input && zip -q -r ../wt_input_1.zip .)
(cd mut_integrate_mut_input && zip -q -r ../mut_input_1.zip .)
rm -rf mut_integrate_wt_input mut_integrate_mut_input
echo "Rebuilt wt_input.zip and mut_input.zip for enzywizard_mut_integrate."
