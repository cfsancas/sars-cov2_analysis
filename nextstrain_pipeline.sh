# Nextstrain SARS-COV-2 pipeline

augur filter \
	--metadata ./data/auspice_metadata.tsv \
	--sequences ./data/sequence.fasta \
	--sequences-per-group 50 \
	--group-by 'lineage' \
	--include-where 'include=yes' \
	--exclude-where 'include=no' \
	--output-metadata ./results/auspice_metadata_filter.tsv \
	--output-sequences ./results/sequence_filter.fasta \
	--output-strains ./results/pass_strains_filter.txt \
	--output-log ./results/filter_log.tsv &&

augur align \
	--sequences ./results/sequence_filter.fasta \
	--output ./results/aligned.fasta \
	--reference-sequence ./data/ref_sequence_ncbi.gb  \
	--remove-reference \
	--nthreads 40 &&

augur tree \
	--alignment results/aligned.fasta \
	--output results/tree_raw.nwk \
	--tree-builder-args '-ninit 10 -n 4' \
	--nthreads 40 &&

## Esta línea se comenta ya que cuando hay secuencias imputadas hay que indicar el clock rate para que no se disparen las fechas
# augur refine --tree results/tree_raw.nwk --alignment results/aligned.fasta --metadata data/auspice_metadata.tsv --output-tree results/tree.nwk --output-node-data results/branch_lengths.json --timetree --coalescent opt --date-confidence --date-inference marginal &&

## Esta es la línea que utilizaba antes de ver los tutoriales de nextstrain apra sars
# augur refine --tree results/tree_raw.nwk --alignment results/aligned.fasta --metadata data/auspice_metadata.tsv --output-tree results/tree.nwk --output-node-data results/branch_lengths.json --timetree --coalescent opt --no-covariance --clock-rate 0.03 &&

augur refine \
	--tree results/tree_raw.nwk \
	--alignment results/aligned.fasta \
	--metadata ./results/auspice_metadata_filter.tsv \
	--output-tree results/tree.nwk \
	--output-node-data results/branch_lengths.json \
	--timetree \
	--divergence-unit mutations \
	--coalescent skyline \
	--date-confidence \
	--clock-rate 0.0008 \
	--clock-std-dev 0.0004 \
	--date-inference marginal \
	--no-covariance &&

# He tenido que quitar --root ya que me daba un error

## esta línea se comenta para evitar un error cuando hay muchas muestras
# augur traits --tree results/tree.nwk --metadata data/auspice_metadata.tsv --output results/traits.json --columns location --confidence &&

# augur frequencies \
# 	--method kde \
# 	--metadata data/auspice_metadata.tsv \
# 	--tree results/tree_raw.nwk \
# 	--min-date 2020-01-01 \
# 	--max-date 2022-05-15 \
# 	--pivot-interval 1 \
# 	--pivot-interval-units weeks \
# 	--narrow-bandwidth 0.05 \
# 	--proportion-wide 0.0 \
# 	--output results/tip-frequencies.json &&

augur ancestral \
	--tree results/tree.nwk \
	--alignment results/aligned.fasta \
	--output-node-data results/nt_muts.json \
	--inference joint \
	--infer-ambiguous &&

# augur distance \ # Da un error al no tener map revisar tutorial
# 	--tree results/test_tree.nwk \
# 	--alignment results/test_aligned.fasta \
# 	--gene-names S \
# 	--compare-to root \
# 	--attribute-name mutations \
# 	--output results/test_distances.json

augur translate \
	--tree results/tree.nwk \
	--ancestral-sequences results/nt_muts.json \
	--reference-sequence ./data/ref_sequence_ncbi.gb \
	--output results/aa_muts.json &&

## Esto es para evitar el error al validad el json: FATAL: Maximum recursion depth reached. You can set the env variable AUGUR_RECURSION_LIMIT to adjust this (current limit: 1000)
export AUGUR_RECURSION_LIMIT=100000

## elimino de los node-data: results/traits.json para que no me de el error: ERROR: meta data file (data/auspice_metadata.tsv) does not exist y he eliminad ese fichero al dar un error al estar hecho con una versión antigua.
augur export v2 \
	--tree results/tree.nwk \
	--metadata ./results/auspice_metadata_filter.tsv \
	--node-data \
		results/branch_lengths.json \
		results/nt_muts.json \
		results/aa_muts.json \
	--lat-longs config/lat_long.tsv \
	--auspice-config config/auspice_config.json \
	--output auspice/SARS-COV-2-$(date -I).json

# results/test_tip-frequencies.json \

## para lanzar auspice
# auspice view --datasetDir auspice

ln -sf auspice/SARS-COV-2-$(date -I).json auspice/SARS-COV-2-latest.json
