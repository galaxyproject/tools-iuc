#!/bin/sh
#
# CHECK-entries.sh
# Marco van Zwetselaar <zwets@kcri.ac.tz>
#
#  Checks that each sequence in the fsa files has an entry in phenotypes.txt,
#  and that every entry in phenotypes has a corresponding sequence in an fsa.
#
#  Writes to stdout the lists of missing entries, as well as close matches
#  based on allele ID or accession.
#

# Echo all identifiers from the phenotypes.txt file, one per line.
phenotype_ids() {
    cut -f1 phenotypes.txt | tail -n +2 | sed -e 's/ *$//'
}

# Echo all sequence identifiers from the *.fsa, one per line.
sequence_ids() {
    fgrep -h '>' *.fsa | sed -e 's/>\([^ ]*\).*/\1/'
}

# Filter stdin for near matches of seqid $1.  A near match is when either
# the allele ID (without accession) or the accession matches.
near_matches() {
    # Horrid REGEX, but we need to escape special chars and split the SEQID
    local REGEX="$(echo "$1" |
        sed -e 's,\([][().+?|*]\),\\\1,g' \
            -e 's,^\([^_]*_[^_]*\)\(_.*\),^\1|\2$,')"
    grep -E "$REGEX" | tr '\n' ' '
}


printf "
===============================================================================
I. Entries in phenotypes.txt with trailing whitespace in their identifier (col 1).
   (This whitespace breaks simple key based lookups.)
-------------------------------------------------------------------------------\n"
cut -f1 phenotypes.txt | tail -n +2 | fgrep ' '


printf "
===============================================================================
II. Entries in phenotypes.txt with no matching sequence in an *.fsa file.
    Second column lists close matches (having identical alleleID or accession).
-------------------------------------------------------------------------------\n"
phenotype_ids | while read SEQID; do
    sequence_ids | fgrep -xq "${SEQID}" ||
        printf "$SEQID\t%s\n" "$(sequence_ids | near_matches "$SEQID")"
done


printf "
===============================================================================
III. Sequences in *.fsa that have no corresponding entry in phenotypes.txt.
     Second column lists close matches (having identical alleleID or accession).
-------------------------------------------------------------------------------\n"
sequence_ids | while read SEQID; do
    phenotype_ids | fgrep -xq "${SEQID}" ||
        printf "$SEQID\t%s\n" "$(phenotype_ids | near_matches "$SEQID")"
done

# vim: sts=4:sw=4:si:et:ai
