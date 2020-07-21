##############
How PuMA Works
##############





The pipeline works in numerous steps. After the sequence data is parsed from the fasta file, each genome is linearlized based on the L1 gene. Next, all main genes are identifed using BLAST (E1, E2, E6, E7, E9, E10, L1, L2). After the main genes are identfied, they are verified by using a database of genes from PaVE by utilizing BLAST and MUSCLE to do sequence identfication, alignment and comparisons. E5 variants (E5_alpha, E5_beta, E5_delta, E5_epsilon, E5_gamma, E5_zeta) are potentially identified by using BLAST to search the genome in the region between the end of E2 and the start of L1. Next, the Upstream Regulatory Region (URR) is identified. Using the URR, E1 and E2 binding sites are identifed via MEME and FIMO. The spliced genes are then identified. The splice acceptor is shared between the E1^E4 and E8^E2 and is embedded with in the E2 gene. MUSCLE is used to align the newly annotated E2 to its closest previously known relative (from BLAST results). For the splice donor site (which is located in E1) a similar approach is used. Once all of the above annotations are created, the various outputs are created. A log file describing the progress of analysis, a 'comma separated values' file (.csv) that contains individual annotation, a 'general features format 3' (.gff3) formated file, a GenBank (.gb) formated file, a PDF file providing a visual representation of the newly annotated genome, and sequin files that are used to streamline the submission process to GenBank are all available after execution. 
