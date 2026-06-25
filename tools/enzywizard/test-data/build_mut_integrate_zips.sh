#!/usr/bin/env bash
set -euo pipefail
base="https://raw.githubusercontent.com/bioinfbrad/enzywizard-mut-integrate/main/examples"
rm -rf wt_input_tmp mut_input_tmp
mkdir -p wt_input_tmp mut_input_tmp

# Download only the reports for 1ZG4_WT and 1ZG6_S70G
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

for f in "${wt_files[@]}"; do
    curl -fsSL "$base/wt_input/$f" -o "wt_input_tmp/$f"
done
for f in "${mut_files[@]}"; do
    curl -fsSL "$base/mut_input/$f" -o "mut_input_tmp/$f"
done

# Remove any extra files (like mut_clean_report or reports for other proteins)
find wt_input_tmp mut_input_tmp -type f -name 'mut_clean_report*.json' -delete

# Create the zip files
rm -f wt_input_1.zip mut_input_1.zip
(cd wt_input_tmp && zip -q -r ../wt_input_1.zip .)
(cd mut_input_tmp && zip -q -r ../mut_input_1.zip .)

rm -rf wt_input_tmp mut_input_tmp
echo "Rebuilt wt_input_1.zip and mut_input_1.zip (only 1ZG4/1ZG6 reports)."