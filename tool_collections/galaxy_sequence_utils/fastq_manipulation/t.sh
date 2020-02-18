TARGET_TOOL_DIRS=("tools/fastq_trimmer_by_quality" "tool_collections/galaxy_sequence_utils/fastq_combiner" "tool_collections/galaxy_sequence_utils/fastq_manipulation" "tool_collections/galaxy_sequence_utils/fastq_groomer" "tool_collections/galaxy_sequence_utils/fastq_filter")

for dir in ${TARGET_TOOL_DIRS[@]}
do
    echo $dir
done
