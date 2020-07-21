######
Output
######


After execution, a folder 'puma_out' will have been created. Within 'puma_out' there will be a folder named 'Short name' (name used in the input file) for each genome in the input file. Within the 'Short name' folder(s) there will be 'for_user' and 'program_files' folders. Also in 'puma-out' is a log file that is by default called 'puma_execution.log'. This log file has potential notes about the annoation of each genome. 


'for_user' contains a .csv file containing indvidual annotations, a .gb file containing annotations in GenBank format, a PDF file that has a visual represntation of the annotated genome, and a 'genbank_submission' folder. 'genbank_submission' has files that will aid in the genbank submission process. 

The 'program_files' folder contains all files PuMA generates and uses during execution.


**For Local Output:**

The 'puma_out' folder will be, by default, in the same folder as all other PuMA files (i.e. 'puma.py' etc.)

**For iMicrobe Output:**

The 'puma_out' folder will be listed under the 'Output' section on the 'Jobs' page. 
