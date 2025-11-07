#!/usr/bin/env bash
set -euo pipefail

# Nextstrain SARS-CoV-2 pipeline processing the full dataset without filtering.

METADATA="./data/auspice_metadata.tsv"
SEQUENCES="./data/sequence.fasta"
REFERENCE="./data/ref_sequence_ncbi.gb"
RESULTS_DIR="./results"
AUSPICE_DIR="./auspice"

mkdir -p "${RESULTS_DIR}" "${AUSPICE_DIR}"

augur align \
        --sequences "${SEQUENCES}" \
        --output "${RESULTS_DIR}/aligned.fasta" \
        --reference-sequence "${REFERENCE}" \
        --remove-reference \
        --nthreads 40 &&

augur tree \
        --alignment "${RESULTS_DIR}/aligned.fasta" \
        --output "${RESULTS_DIR}/tree_raw.nwk" \
        --tree-builder-args '-ninit 10 -n 4' \
        --nthreads 40 &&

# Refine the tree using the complete metadata without applying any filtering.
augur refine \
        --tree "${RESULTS_DIR}/tree_raw.nwk" \
        --alignment "${RESULTS_DIR}/aligned.fasta" \
        --metadata "${METADATA}" \
        --output-tree "${RESULTS_DIR}/tree.nwk" \
        --output-node-data "${RESULTS_DIR}/branch_lengths.json" \
        --timetree \
        --divergence-unit mutations \
        --coalescent skyline \
        --date-confidence \
        --clock-rate 0.0008 \
        --clock-std-dev 0.0004 \
        --date-inference marginal \
        --no-covariance &&

augur ancestral \
        --tree "${RESULTS_DIR}/tree.nwk" \
        --alignment "${RESULTS_DIR}/aligned.fasta" \
        --output-node-data "${RESULTS_DIR}/nt_muts.json" \
        --inference joint \
        --infer-ambiguous &&

augur translate \
        --tree "${RESULTS_DIR}/tree.nwk" \
        --ancestral-sequences "${RESULTS_DIR}/nt_muts.json" \
        --reference-sequence "${REFERENCE}" \
        --output "${RESULTS_DIR}/aa_muts.json" &&

# Increase recursion limit to prevent Augur export errors in large datasets.
export AUGUR_RECURSION_LIMIT=100000

augur export v2 \
        --tree "${RESULTS_DIR}/tree.nwk" \
        --metadata "${METADATA}" \
        --node-data \
                "${RESULTS_DIR}/branch_lengths.json" \
                "${RESULTS_DIR}/nt_muts.json" \
                "${RESULTS_DIR}/aa_muts.json" \
        --lat-longs config/lat_long.tsv \
        --auspice-config config/auspice_config.json \
        --output "${AUSPICE_DIR}/SARS-COV-2-$(date -I).json"

ln -sf "SARS-COV-2-$(date -I).json" "${AUSPICE_DIR}/SARS-COV-2-latest.json"
