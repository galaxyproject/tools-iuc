#!/bin/bash -e

# Script called by bbnorm.xml to deterministically pre-estimate the amount of memory to be 
# used by the tool.

estimate_required_memory () {
    loglog_file=$1
    num_hashes=$2
    bits_per_kmer=$3
    memory_filled_pct=$4

    cardinality=$(cat $loglog_file | grep Cardinality | awk '{print $2}')

    # Example memory calculation:
    # 3hashes * 16bits/kmer * 1000000000kmers * 100/50 pct_load / 8000000bits/MB = 12 GB to achieve a table under 50% full + 256 MB extra buffer
    required_RAM_MB=$(( ${num_hashes} * ${bits_per_kmer} * $cardinality * 100 / ${memory_filled_pct} / 8000000 + 256))

    echo Estimated kmer cardinality of the input: $cardinality
    echo Estimated RAM requirement: ${required_RAM_MB}MB
}


check_memory_limits () {
    # Check if an external constraint limits the memory, and if so, whether we are below that limit.
    # If not, throw an out of memory error.
    if [ $required_RAM_MB -gt "${GALAXY_MEMORY_MB:-9192}" ]; then
        >&2 echo Your Galaxy configuration limits RAM usage to "${GALAXY_MEMORY_MB:-9192}"MB, whereas you need at least ${required_RAM_MB}MB to accurately perform this task. Increase \"Hashing parameters\" ">" \"Percent memory filled\" or contact your Galaxy administrator.
        exit 42
    fi
    
    if [ -n "${_JAVA_OPTIONS}" ]; then
        source calcmem.sh
        parseXmx "${_JAVA_OPTIONS}"
        echo $z
	    if [[ $set == 1 ]]; then
	        if [[ ${z:0-1} == 'm' ]]; then
	            JAVA_OPT_MEMORY_MB=${z:4:-1}
	        elif [[ ${z:0-1} == 'g' ]]; then
	            JAVA_OPT_MEMORY_MB=$((${z:4:-1}*1024))
	        fi
	        
            if [ $required_RAM_MB -gt "$JAVA_OPT_MEMORY_MB" ]; then
                >&2 echo Your Java configuration limits RAM usage to "$JAVA_OPT_MEMORY_MB"MB, whereas you need at least ${required_RAM_MB}MB to accurately perform this task. Increase \"Hashing parameters\" ">" \"Percent memory filled\" or contact your Galaxy administrator.
                exit 42
            fi
        fi
    fi
    export _JAVA_OPTIONS="$_JAVA_OPTIONS -Xmx${required_RAM_MB}m -Xms${required_RAM_MB}m"
}


