# sequence_handling
create, update and filtering fastq database

### Operating system: Linux 

### External Tools:  

    - [Cutadapt](https://cutadapt.readthedocs.io/en/stable/): for trimming 

    - [NGmerge](https://github.com/jsh58/NGmerge): for merging 

    - [Seqtk](https://github.com/lh3/seqtk): for filtering based on min_amplicon_length 

### Database: [sqlite3](https://www.sqlite.org/index.html) 


### Installation: 
1. Clone the source code from github (with git clone) 
2. Install cutadapt and seqtk: 
``` 
    sudo apt install cutadapt 
    sudo apt install seqtk 

```
3. NGmerge is already provided in the project in “Tools” folder. Copy the NGmerge application to another directory, in case of permission error. If it still produces some other error, manuel Install NGmerge (see github)
** important: adjust path to NGmerge in utils.py**  
     

4. Download sqlite3 version 3.41.2  
     
5. Create python env is necessary: ``` python -m venv <env_name> ``` 

6. Activate env: ```  source <env_name>/bin/activate ``` 

7. Go to project directory 

8. Install python packages: ``` pip install –r requirements.txt  ``` 
     
Test Data: [Onedrive](https://pages.github.com/](https://1drv.ms/f/c/33A4E3C9AE76E077/Erj5omHYGUhIsYYfL1u20n4Bp7eXtejy3jaN2rUlngI8Tg?e=Xa3VAa).
Gained from: https://www.applied-maths.com/download/fastq-files  
     
## Running Test Data 

1. Download Test Data, unzipped data 
2. Got to the project directory: ``` cd <dir>/sequence_handling ``` 
3. Adjust parameter/directories in test.py (Test.py will call main.py)  
4. Activate python env 
5. Run test.py -> ``` python test.py ``` 

Run time= O(merged_fastq_reads_len * db_length) 

Test Data= around 30-40 minutes 
