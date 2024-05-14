
from main import *

'''
Test scenario: 
- creating db with fastq files in db input folder (dbFolder)
- filtering with input_fastq1 and input_fastq2 

Args:
    input_fastq1 (str): Path to the first input FASTQ file.
    input_fastq2 (str): Path to the second input FASTQ file.
    threshold (int): Threshold value for filtering sequences.
    output_fasta (str): Path to the output FASTA file.
    outdir (str): Directory to store output files.
    db_name (str): Name of the database.
    db_folder (str, optional): Path to the database folder. Defaults to None.
    min_overlap(int, default:200): Minimum overlap length 
    min_amplicon_length(int, default:400): Minimum Amplicon length
    fixed_sequence1 (str, optional): Fixed sequence for input FASTQ 1. Defaults to None.
    fixed_sequence2 (str, optional): Fixed sequence for input FASTQ 2. Defaults to None.
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
