############################
Local Usage (MacOS or Linux)
############################



Due to dependencies, PuMA only runs on MacOS or Linux operating systems. Refer to the 'Dependency Information' section for more information. 

After the dependencies have been installed, from `Github <https://github.com/KVD-lab/puma>`_ download the 'data_dir' folder and within the 'scripts' directory download 'run_puma.py' and 'puma.py'.

To execute PuMA, make sure 'run_puma.py' is excecutable:
::

    chmod +x run_puma.py

Minimum PuMA execuation command:
::
    
	./run_puma.py -i HPV16REF.fa -D warning

All files and folders need to be in the same location or absolute path needs to be used. For example the input file above (HPV16REF) and the 'data_dir' folder need to be in the same folder as 'run_puma.py' and 'puma.py'. 



You can always use '-h' or '--help' on 'run_puma.py' to get a list and description of all inputs. 

**After execution, it is encouraged to either move, rename, or delete the output folder 'puma_out' before running PuMA again. This way, there are no issues with the existing 'puma_out' folder and the new one that will be generated with each execution.** 
