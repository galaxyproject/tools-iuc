set.seed(910401)

uniq <- sample(c('A', 'T', 'C', 'G'), 7000, replace = T)
dupl <- sample(c('A', 'T', 'C', 'G'), 500, replace = T)
# duplicaiton will be imperfect
mutated_dupl <- dupl
mutated_dupl[round(runif(20, 1, 500) )] <- sample(c('A', 'T', 'C', 'G'), 20, replace = T)

m_haplotype <- c(uniq, dupl, mutated_dupl, mutated_dupl[1:50])
p_haplotype <- c(uniq, dupl, mutated_dupl, mutated_dupl[1:50])

# create some heterozygous sites
p_haplotype[round(runif(100, 1, 8000) )] <- sample(c('A', 'T', 'C', 'G'), 100, replace = T)
# mean(p_haplotype == m_haplotype)
# 99% same haplotypes

m_sequencing <- rep(m_haplotype, 25)
p_sequencing <- rep(p_haplotype, 25)

# create some sequencing errors
m_sequencing[round(runif(200, 1, 200000) )] <- sample(c('A', 'T', 'C', 'G'), 200, replace = T)
p_sequencing[round(runif(200, 1, 200000) )] <- sample(c('A', 'T', 'C', 'G'), 200, replace = T)

# mean(m_sequencing == p_sequencing)
# 98.9% same -> 1% het, 0.1% seq errors

m_header <- '>maternal_hapl_sequences'
m_fasta_entry <- paste0(m_sequencing, collapse = '')
p_header <- '>paternal_hapl_sequences'
p_fasta_entry <- paste0(p_sequencing, collapse = '')

writeLines(c(m_header, m_fasta_entry, p_header, p_fasta_entry), 'smudgeplot_test_data.fasta')