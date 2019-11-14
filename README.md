# PuMA Release 1.1 

**PuMA is now available to access through [iMicrobe](https://www.imicrobe.us).**

Papillomavirus genome annotation tool

Release 1.1 (11/14/19) added .fsa & .tbl formatted output for help with genbank submission process, significantly more robust overall
# Authors

Josh Pace, Ken Younes-Clark, Cordell Freeman, Koenraad Van Doorslaer 

University of Arizona, KVD Lab & Hurwitz Lab

# Formatting Input FASTA File
    
    >Short name|Full Name
    Sequence


Short name is the abbreviation or accession number you want for output files (e.g. HPV16)

Full name is what will be printed to the screen (e.g. Human papillomavirus 16)

# PuMA Output Files

Within the puma_out folder, there are two folders, for_user and program_files. Within the for_user folder there are:
* csv file containing start, stop, nucleotide and protein sequences
* log file that has information from execution
* a graphical output of the open reading frames and miscellaneous features 
* genbank file 
* gff3 file 
* .fsa file
* .tbl file

# Dependencies 

Please install the following:

* Python 3.x
* Biopython (pip)
* NCBI BLAST+ 2.7.x (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* MEME, FIMO (https://meme-suite.org/)
* MUSCLE (https://www.drive5.com/muscle)
* pandas (pip)
* matplotlib (pip)

# PuMA Poly Release 0.1 

Polyomavirus genome annotation tool 
Annoates VP1, VP2, small t and Large T proteins

# Authors

Josh Pace, Ken Younes-Clark, Koenraad Van Doorslaer 

University of Arizona, KVD Lab & Hurwitz Lab

Release 0.1 (7/25/19) output is printed genome information and csv

# Dependencies 

Please install the following:

* Python 3.x
* Biopython (pip)
* NCBI BLAST+ 2.7.x (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* MUSCLE (https://www.drive5.com/muscle)

# Formatting Input FASTA File
    
    >Short name|Full Name
    Sequence


Short name is the abbreviation or accession number you want for output files (e.g. SV40a)

Full name is what will be printed to the screen (e.g. simian virus 40)
