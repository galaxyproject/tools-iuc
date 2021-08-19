#!/bin/sh

NEXTCLADE_INPUTS_DIR=db
NEXTCLADE_INPUTS_URL_BASE="https://raw.githubusercontent.com/nextstrain/nextclade/master/data/sars-cov-2"
echo "[ INFO] ${0}:${LINENO}: Downloading Nextclade input data:"
echo "[ INFO] ${0}:${LINENO}:   From: \"${NEXTCLADE_INPUTS_URL_BASE}\""
echo "[ INFO] ${0}:${LINENO}:   To:   \"${NEXTCLADE_INPUTS_DIR}\""
mkdir -p "${NEXTCLADE_INPUTS_DIR}"
pushd "${NEXTCLADE_INPUTS_DIR}" >/dev/null
  curl -fsSLOJ --write-out "[ INFO] curl: %{url_effective}\n" "${NEXTCLADE_INPUTS_URL_BASE}/reference.fasta"
  curl -fsSLOJ --write-out "[ INFO] curl: %{url_effective}\n" "${NEXTCLADE_INPUTS_URL_BASE}/genemap.gff"
  curl -fsSLOJ --write-out "[ INFO] curl: %{url_effective}\n" "${NEXTCLADE_INPUTS_URL_BASE}/tree.json"
  curl -fsSLOJ --write-out "[ INFO] curl: %{url_effective}\n" "${NEXTCLADE_INPUTS_URL_BASE}/qc.json"
  curl -fsSLOJ --write-out "[ INFO] curl: %{url_effective}\n" "${NEXTCLADE_INPUTS_URL_BASE}/primers.csv"
popd >/dev/null
