# Check an OTU table before mumu reads it, and write a clean copy to stdout.
#
# mumu reads every column after the first as a whole number and does not warn
# when it cannot: 1e+05 is read as 1, 0.9 as 0 and text as 0. This script:
#   - removes Windows line endings and skips empty lines;
#   - rewrites whole numbers in scientific notation (1e+05) or with a zero
#     decimal part (12.0) as plain whole numbers;
#   - stops with an error for decimals, negative numbers, empty cells, text and
#     lines with a different number of columns from the first line.
# Use: awk -f check_otu_table.awk otu_table.tsv > clean_table.tsv

BEGIN {
    FS = "\t"
    OFS = "\t"
    converted = 0
}

function fail(message) {
    printf "Error in the OTU or ASV table, line %d: %s\n", FNR, message > "/dev/stderr"
    failed = 1
    exit 1
}

function column_name(i) {
    return "column " i " (" header[i] ")"
}

{ sub(/\r$/, "") }

/^[ \t]*$/ { next }

!seen_header {
    seen_header = 1
    columns = NF
    if (columns < 2) {
        fail("the first line has only one column. The table must be tab-separated, with OTU or ASV names in the first column and one column per sample.")
    }
    for (i = 1; i <= NF; i++) {
        header[i] = $i
    }
    print
    next
}

{
    if (NF != columns) {
        fail(sprintf("this line has %d columns, but the first line has %d. Check that the file is tab-separated and that every line has a value for every sample.", NF, columns))
    }
    for (i = 2; i <= NF; i++) {
        value = $i
        gsub(/^ +| +$/, "", value)
        if (value ~ /^[0-9]+$/) {
            $i = value
        } else if (value ~ /^[0-9]+\.0*$/) {
            sub(/\..*/, "", value)
            $i = value
            converted++
        } else if (value ~ /^[+]?([0-9]+\.?[0-9]*|\.[0-9]+)[eE][+-]?[0-9]+$/) {
            number = value + 0
            if (number != int(number)) {
                fail(sprintf("%s holds %s, which is not a whole number. mumu needs read counts, not proportions or decimals.", column_name(i), value))
            }
            if (number > 9007199254740992) {
                fail(sprintf("%s holds %s, which is too large to read exactly.", column_name(i), value))
            }
            $i = sprintf("%.0f", number)
            converted++
        } else if (value ~ /^-/) {
            fail(sprintf("%s holds %s. Read counts cannot be negative.", column_name(i), value))
        } else if (value ~ /^[+]?([0-9]+\.[0-9]*|\.[0-9]+)$/) {
            fail(sprintf("%s holds %s, which is not a whole number. mumu needs read counts, not proportions or decimals.", column_name(i), value))
        } else if (value == "") {
            fail(sprintf("%s is empty. Every sample needs a read count (0 if absent).", column_name(i)))
        } else {
            fail(sprintf("%s holds the text '%s'. All columns after the first must hold read counts, so remove sequence or taxonomy columns first (for example with the Cut tool).", column_name(i), substr(value, 1, 40)))
        }
    }
    rows++
    print
}

END {
    if (failed) {
        exit 1
    }
    if (!seen_header) {
        printf "Error in the OTU or ASV table: the file is empty.\n" > "/dev/stderr"
        exit 1
    }
    if (rows == 0) {
        printf "Error in the OTU or ASV table: there is a header line but no OTUs or ASVs.\n" > "/dev/stderr"
        exit 1
    }
    if (converted > 0) {
        printf "OTU or ASV table: %d values written in scientific notation or with a decimal point (for example 1e+05 or 12.0) were changed to whole numbers.\n", converted > "/dev/stderr"
    }
}
