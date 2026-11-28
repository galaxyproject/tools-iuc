# Convert FasTAN annotation (.1ano, read through ONEview) to BED6 in scaffold coordinates.
#
# .1ano schema (FasTAN ANO.c, anoSchemaText):
#     O S 1 6 STRING              id for a scaffold
#     D G 1 3 INT                 gap of given length
#     D C 1 3 INT                 contig of given length
#     O M 3 3 INT 3 INT 3 INT     scaffold index, beg,end pair
#     D L 1 6 STRING              optional label for preceding M
#     D X 1 3 INT                 optional score for preceding M
#
# M holds a 0-based scaffold index (the order of the S lines) and scaffold-relative coordinates.
# This script reads them directly. ANOtoBED is not used, because it reports contig-relative
# coordinates under the scaffold name, which is wrong after an internal N gap.
#
# ONEcode ASCII writes each STRING field with a length prefix, so the value of S and L is $3.
# INT fields have no prefix, so the value of X is $2.
#
# L holds the repeat unit length. X is the divergence (about a percentage; 0 is a perfect array).
# Column 5 is 1000 - 10 * divergence, clamped to 0-1000, so a higher score means a purer array.
# A record is written on the next M or at END, because optional records (such as P) can follow M.
BEGIN { OFS = "\t"; nscaf = -1; have = 0 }

function flush() {
    if (!have) return
    d = score + 0; if (d < 0) d = 0; if (d > 100) d = 100
    print name[s], b, e, "u" unit, int(1000 - 10 * d), "."
    have = 0
}

$1 == "S" { nscaf++; name[nscaf] = $3; next }
$1 == "M" { flush(); s = $2 + 0; b = $3; e = $4; unit = "?"; score = 0; have = 1; next }
$1 == "L" && have { unit = $3; next }
$1 == "X" && have { score = $2; next }
END { flush() }
