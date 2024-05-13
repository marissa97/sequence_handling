
from main import *

'''
Test scenario: 
- creating db with fastq files in db input folder (dbFolder)
- filtering with input_fastq1 and input_fastq2 
'''

dbFolder= "/home/marissa/Documents/dkms/input_DB"
outdir="/home/marissa/Documents/dkms/output"
input_fastq1="/home/marissa/Documents/dkms/input_filtering/ERR101900_1.fastq.gz"
input_fastq2="/home/marissa/Documents/dkms/input_filtering/ERR101900_2.fastq.gz"
db_name="test.db"
threshold= 10 ## levenshtein distance threshold (value to 10, for test faster)
min_overlap=15
min_amplicon_length=100


#conn, cursor= activateDB(db_name)
#mainCreateDB(conn, cursor, inputFolder, 15, 100,outdir, threshold=3)

main(input_fastq1, input_fastq2,  outdir, db_name,  min_overlap, min_amplicon_length, threshold, dbFolder)
