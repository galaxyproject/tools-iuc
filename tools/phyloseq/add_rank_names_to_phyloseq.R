#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("phyloseq"))
suppressPackageStartupMessages(library("tidyverse"))

# Option parsing
option_list <- list(
    make_option(c("--input"),
        action = "store", dest = "input",
        help = "Input file containing a phyloseq object"
    ),
    make_option(c("--output"), action = "store", dest = "output", help = "Output file for the updated phyloseq object")
)

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser, positional_arguments = TRUE)
opt <- args$options

cat("Input file: ", opt$input, "\n")
cat("Output file: ", opt$output, "\n")

# Lade das Phyloseq-Objekt
physeq <- readRDS(opt$input)

# Überprüfen, ob das Phyloseq-Objekt erfolgreich geladen wurde
if (is.null(physeq)) {
    stop("Error: Failed to load the Phyloseq object. Check the input file.")
}

cat("Phyloseq object successfully loaded.\n")
cat("Class of loaded object: ", class(physeq), "\n")

# Überprüfe das aktuelle Taxonomy Table
cat("Current tax_table:\n")
print(tax_table(physeq))

# Zuordnung der Rangnamen
rank_names <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") # Anpassen je nach Bedarf

# Füge eine leere Spalte für Species hinzu, falls sie fehlt
if (ncol(tax_table(physeq)) == 6) {
    tax_table(physeq) <- cbind(tax_table(physeq), Species = NA)
}

# Überprüfen, ob die Anzahl der Spalten mit der Anzahl der Rangnamen übereinstimmt
if (ncol(tax_table(physeq)) != length(rank_names)) {
    stop("Error: Number of columns in tax_table does not match the length of rank_names.")
}

# Setzen der Spaltennamen
colnames(tax_table(physeq)) <- rank_names

# Bestätige die Änderungen
cat("Updated tax_table:\n")
print(tax_table(physeq))

# Extrahiere das erste Zeichen aus dem ersten Eintrag des tax_table (z.B. Kingdom)
first_char <- substr(tax_table(physeq)[1, 1], 1, 1)
cat("Extracted first character: ", first_char, "\n")

# Speichere das aktualisierte Phyloseq-Objekt
saveRDS(physeq, file = opt$output, compress = TRUE)
cat("Updated Phyloseq object saved to: ", opt$output, "\n")

# Gib das erste Zeichen zurück (es wird später in der XML-Datei verwendet)
cat(first_char, "\n")
