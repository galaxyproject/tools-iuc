# Rewrite selected fields of every @RG line of a SAM header.
#
# Usage: awk -f edit_read_groups.awk edits.tsv - < header.sam
#
# edits.tsv holds one "TAG<TAB>VALUE" record per field to set. The values are
# only ever read as data and are inserted by plain string concatenation, so no
# user-supplied string can reach a shell, be evaluated as an awk expression, or
# be interpreted as a sub()/gsub() replacement.
#
# Every @RG record is edited in place: existing TAG:VALUE fields are replaced
# and missing ones are appended. ID (and any other field not listed in
# edits.tsv) is left untouched, which is what keeps per-lane read groups
# distinct while SM/LB become per-sample constants.
BEGIN {
    FS = "\t"
    OFS = "\t"

    # Read the edits up front rather than as awk's first input file: the usual
    # NR == FNR idiom silently treats the start of the header as edits when the
    # table is empty.
    edits_file = ARGV[1]
    ARGV[1] = ""
    while ((getline line < edits_file) > 0) {
        if (split(line, field, "\t") < 2 || field[1] == "") {
            continue
        }
        if (field[2] == "") {
            print "edit_read_groups: no value given for @RG field " field[1] > "/dev/stderr"
            failed = 1
            exit 1
        }
        if (!(field[1] in edit)) {
            order[++n_edits] = field[1]
        }
        edit[field[1]] = field[2]
    }
    close(edits_file)

    if (n_edits == 0) {
        print "edit_read_groups: no @RG fields to set" > "/dev/stderr"
        failed = 1
        exit 1
    }
}

$1 == "@RG" {
    n_rg++
    split("", seen)
    for (i = 2; i <= NF; i++) {
        tag = substr($i, 1, 2)
        if (substr($i, 3, 1) == ":" && tag in edit) {
            $i = tag ":" edit[tag]
            seen[tag] = 1
        }
    }
    for (i = 1; i <= n_edits; i++) {
        if (!(order[i] in seen)) {
            $(NF + 1) = order[i] ":" edit[order[i]]
        }
    }
    print
    next
}

{
    print
}

END {
    if (failed) {
        exit 1
    }
    # Silently returning the header unchanged is the worst outcome here: the job
    # goes green and the user believes the sample was renamed.
    if (n_rg == 0) {
        print "edit_read_groups: the input header contains no @RG records, so there is nothing to edit" > "/dev/stderr"
        exit 1
    }
}
