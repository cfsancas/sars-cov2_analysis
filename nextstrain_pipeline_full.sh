#!/usr/bin/env bash
set -euo pipefail

# Nextstrain SARS-CoV-2 pipeline processing the full dataset without filtering.

METADATA="./data/auspice_metadata.tsv"
SEQUENCES="./data/sequence.fasta"
REFERENCE="./data/ref_sequence_ncbi.gb"
RESULTS_DIR="./results_full"
AUSPICE_DIR="./auspice_full"

mkdir -p "${RESULTS_DIR}" "${AUSPICE_DIR}"

time augur filter \
  --sequences "${SEQUENCES}" \
  --metadata "${METADATA}" \
  --output-sequences "${RESULTS_DIR}/filtered.fasta" \
  --output-metadata  "${RESULTS_DIR}/metadata.filtered.tsv" \
  --min-length 29500 \
  --max-ambiguous 3000 \
  --exclude-ambiguous-dates-by any \
  --min-date 2019-12-01 \
  --prune-duplicates \
  --subsample-max-sequences 0

time augur mask \
  --sequences "${RESULTS_DIR}/filtered.fasta" \
  --output "${RESULTS_DIR}/masked.fasta" \
  --mask-from-beginning 100 \
  --mask-from-end 100

time augur align \
    --sequences "${RESULTS_DIR}/masked.fasta" \
    --output "${RESULTS_DIR}/aligned.fasta" \
    --reference-sequence "${REFERENCE}" \
    --remove-reference \
    --method "nextalign" \
    --nthreads 40

time augur tree \
    --alignment "${RESULTS_DIR}/aligned.fasta" \
    --output "${RESULTS_DIR}/tree_raw.nwk" \
    --method iqtree \
    --tree-builder-args '-seed 12345  -m GTR+G' \
    --nthreads 40 &&

# Refine the tree using the complete metadata without applying any filtering.
time augur refine \
    --tree "${RESULTS_DIR}/tree_raw.nwk" \
    --alignment "${RESULTS_DIR}/aligned.fasta" \
    --metadata "${RESULTS_DIR}/metadata.filtered.tsv" \
    --output-tree "${RESULTS_DIR}/tree.nwk" \
    --output-node-data "${RESULTS_DIR}/branch_lengths.json" \
    --timetree \
    --divergence-unit mutations \
    --coalescent skyline \
    --date-confidence \
    --clock-rate 0.0008 \
    --clock-std-dev 0.0004 \
    --clock-filter-iqd 4 \
    --date-inference marginal \
    --no-covariance

time augur ancestral \
    --tree "${RESULTS_DIR}/tree.nwk" \
    --alignment "${RESULTS_DIR}/aligned.fasta" \
    --output-node-data "${RESULTS_DIR}/nt_muts.json" \
    --inference joint \
    --infer-ambiguous

time augur translate \
    --tree "${RESULTS_DIR}/tree.nwk" \
    --ancestral-sequences "${RESULTS_DIR}/nt_muts.json" \
    --reference-sequence "${REFERENCE}" \
    --genes S ORF1a ORF1b N E M \
    --output "${RESULTS_DIR}/aa_muts.json"

time augur traits \
  --tree "${RESULTS_DIR}/tree.nwk" \
  --metadata "${RESULTS_DIR}/metadata.filtered.tsv" \
  --columns region country division \
  --output-node-data "${RESULTS_DIR}/traits.json" \
  --confidence

  time augur frequencies \
    --tree "${RESULTS_DIR}/tree.nwk" \
    --metadata "${RESULTS_DIR}/metadata.filtered.tsv" \
    --output "${RESULTS_DIR}/frequencies.json"

# Increase recursion limit to prevent Augur export errors in large datasets.
export AUGUR_RECURSION_LIMIT=100000

augur export v2 \
        --tree "${RESULTS_DIR}/tree.nwk" \
        --metadata "${RESULTS_DIR}/metadata.filtered.tsv" \
        --node-data \
                "${RESULTS_DIR}/branch_lengths.json" \
                "${RESULTS_DIR}/nt_muts.json" \
                "${RESULTS_DIR}/aa_muts.json" \
                "${RESULTS_DIR}/traits.json" \
                "${RESULTS_DIR}/frequencies.json" \
        --lat-longs config/lat_long.tsv \
        --geo-resolutions country division location \
        --auspice-config config/auspice_config.json \
        --output "${AUSPICE_DIR}/SARS-COV-2-$(date -I).json"

ln -sf "SARS-COV-2-$(date -I).json" "${AUSPICE_DIR}/SARS-COV-2-latest.json"
