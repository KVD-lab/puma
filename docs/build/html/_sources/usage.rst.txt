#####
Usage
#####

Make sure that `Biopython <https://biopython.org>`_, `BLAST <https://blast.ncbi.nlm.nih.gov/    Blast.cgi?PAGE_TYPE=BlastDocs>`_, `MEME, FIMO <https://meme-suite.org/>`_, `MUSCLE <https://www.drive5.com/muscle>`_, pandas (pip),  matplotlib (pip) and Python 3.x are installed.

From `Github <https://github.com/KVD-lab/puma>`_ download the 'data_dir' directory and within the 'scripts' directory download 'run_puma.py' and 'puma.py'.

To execute PuMA, make sure 'run_puma.py' is excecutable:
::

    chmod +x run_puma.py

PuMA execuation command:
::
    
	./run_puma.py -i HPV16REF.fa  -d data_dir -D warning

All files and folders need to be in the same location or absolute path needs to be used. For example the input file above (HPV16REF) and the 'data_dir' need to be in the same folder as 'run_puma.py' and 'puma.py'. 
