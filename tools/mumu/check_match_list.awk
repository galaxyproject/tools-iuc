# Check a match list before mumu reads it, and write a clean copy to stdout.
#
# This script:
#   - removes Windows line endings and skips empty lines;
#   - skips a header line (a first line whose third column is not a number);
#   - stops with an error for lines that do not have three columns, for
#     similarities that are not numbers from 0 to 100, and for match lists
#     where all similarities are 1 or less (fractions instead of percentages,
#     which mumu would accept and then merge nothing).
# Use: awk -f check_match_list.awk match_list.tsv > clean_match_list.tsv

BEGIN {
    FS = "\t"
    OFS = "\t"
    number_pattern = "^[+]?([0-9]+\\.?[0-9]*|\\.[0-9]+)([eE][+-]?[0-9]+)?$"
}

function fail(message) {
    printf "Error in the match list, line %d: %s\n", FNR, message > "/dev/stderr"
    failed = 1
    exit 1
}

{ sub(/\r$/, "") }

/^[ \t]*$/ { next }

{
    if (NF != 3) {
        fail(sprintf("this line has %d columns. The match list needs exactly three: query OTU or ASV name, hit OTU or ASV name, percent similarity. For BLASTN standard tabular output, keep columns 1 to 3 with the Cut tool.", NF))
    }
    value = $3
    gsub(/^ +| +$/, "", value)
    if (value !~ number_pattern) {
        if (lines == 0 && !skipped_header) {
            skipped_header = 1
            next
        }
        fail(sprintf("the similarity in column 3 is '%s', which is not a number.", substr(value, 1, 40)))
    }
    similarity = value + 0
    if (similarity > 100) {
        fail(sprintf("the similarity in column 3 is %s, which is more than 100. Column 3 must be the percent similarity.", value))
    }
    if (similarity > maximum) {
        maximum = similarity
    }
    $3 = value
    lines++
    print
}

END {
    if (failed) {
        exit 1
    }
    if (lines == 0) {
        printf "Error in the match list: there are no matches.\n" > "/dev/stderr"
        exit 1
    }
    if (maximum <= 1) {
        printf "Error in the match list: all similarities are 1 or less. They look like fractions (0.97) instead of percentages (97). For VSEARCH, use the id user field, which is a percentage.\n" > "/dev/stderr"
        exit 1
    }
    if (skipped_header) {
        printf "Match list: the first line was taken as a header and skipped.\n" > "/dev/stderr"
    }
}
