The paired FASTQ fixtures simulate error-free Illumina reads from the
`rep14a_2_ORF1(pKQ10)_EFU01917` sequence in `test-db/Rep_trans.fsa`.
Each pair has 150 bp reads from opposite ends of a 300 bp fragment. Fragment
starts cover all positions from 0 to 102, giving 103 pairs and complete coverage
of the 402 bp reference when both mates are used. Quality scores are Phred 40.

Together the mates produce one hit with 100% identity and coverage. Neither
mate alone passes the default identity threshold. This makes the tests detect
an accidentally omitted mate. The plain and gzip files contain the same reads.
