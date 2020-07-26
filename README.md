[![Documentation Status](https://readthedocs.org/projects/puma-docs/badge/?version=latest)](https://puma-docs.readthedocs.io/en/latest/?badge=latest) (click to be directed to the Read the Docs documentation (more detailed)) 


# PuMA Release 1.2.0 

**PuMA is now available to access through [iMicrobe](https://www.imicrobe.us).**

**PuMA is now available to access through [Docker](https://hub.docker.com/r/kvdlab/puma)**

Papillomavirus genome annotation tool. 

Release 1.2.0 (7/24/2020) added gene alignment verification, ability to handle a multi genome input fasta file, and updated databases. 
# Authors

Josh Pace, Ken Younes-Clark, Cordell Freeman, Koenraad Van Doorslaer 

University of Arizona, KVD Lab & Hurwitz Lab

# Formatting Input FASTA File
    
    >Short name|Full Name
    Sequence


Short name is the abbreviation or accession number you want for output files (e.g. HPV16)

Full name is what will be printed to the screen (e.g. Human papillomavirus 16)

For a multi genome file, follow the same naming convention as above. 

# PuMA Output Files

Within the puma_out folder, there is a log file that has information from execution and there will be a folder for each genome named with the "Short name". Within each "Short name" folder there are 'for_user' and 'program_files' folders. Within the for_user folder there are:
* csv file containing start, stop, nucleotide and protein sequences
* a graphical output of the open reading frames and miscellaneous features 
* genbank file 
* gff3 file 
* .fsa file
* .tbl file

# Dependencies 

Please install the following:

* [Python 3.x](https://www.python.org/downloads/)
* [Biopython (pip)](https://biopython.org/wiki/Download)
* [NCBI BLAST+ 2.7.x](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [MUSCLE](https://www.drive5.com/muscle)
* [MEME, FIMO](http://meme-suite.org/doc/install.html?man_type=web)
* [matplotlib (pip)](https://matplotlib.org/users/installing.html)
* [pandas (pip)](https://pandas.pydata.org)

# PuMA Poly Release 0.1 

Polyomavirus genome annotation tool 
Annoates VP1, VP2, small t and Large T proteins

# Authors

Josh Pace, Ken Younes-Clark, Koenraad Van Doorslaer 

University of Arizona, KVD Lab & Hurwitz Lab

Release 0.1 (7/25/19) output is printed genome information and csv

# Dependencies 

Please install the following:

* [Python 3.x](https://www.python.org/downloads/)
* [Biopython (pip)](https://biopython.org/wiki/Download)
* [NCBI BLAST+ 2.7.x](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [MUSCLE](https://www.drive5.com/muscle/downloads.htm)

# Formatting Input FASTA File
    
    >Short name|Full Name
    Sequence


Short name is the abbreviation or accession number you want for output files (e.g. SV40a)

Full name is what will be printed to the screen (e.g. simian virus 40)
