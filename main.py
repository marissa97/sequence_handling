import os
import sys
fpath = os.path.join(os.path.dirname(__file__), 'code')
sys.path.append(fpath)
print(sys.path)

from task1 import *
from task2 import *

# main function
def main( input_fastq1, input_fastq2, outdir, db_name,  min_overlap=200, min_amplicon_length=400,threshold=3, db_folder=None, fixed_sequence1=None, fixed_sequence2=None):

    """
        Main function to process input FASTQ files, filter sequences, and update database.

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
    """
    try:
        # Create database if db_folder is provided
        conn, cursor= activateDB(db_name)
        if db_folder is not None:
            mainCreateDB(conn, cursor, db_folder,  min_overlap, min_amplicon_length,outdir, threshold,fixed_sequence1, fixed_sequence2)

        # Check if input FASTQ files belong to the same sample
        input_fastq1_name = os.path.basename(input_fastq1).split("/")[len(os.path.basename(input_fastq1).split("/"))-1].split("_")[0]
        input_fastq2_name = os.path.basename(input_fastq2).split("/")[len(os.path.basename(input_fastq2).split("/"))-1].split("_")[0]

        if input_fastq1_name != input_fastq2_name:
            print("samples name are different, stopping... ")
            return
    
        results= checkExistFastqFile(cursor, input_fastq1_name)
        
        if len(results) >0: 
            print("the fastq file is already found in the database, stopping ...")
            return

        
        newoutdir= outdir +"/"+ input_fastq1_name
        if not os.path.exists(newoutdir):
            # Create the directory
            os.mkdir(newoutdir)
    

        filter_sequences(conn,cursor,input_fastq1_name,input_fastq1, input_fastq2, threshold, newoutdir, min_overlap, min_amplicon_length,fixed_sequence1, fixed_sequence2 )
        
    except Exception as e:
        print(f"Error occurred: {e}")

    finally:
        # Close database connection
        if 'conn' in locals():
            conn.close()


if __name__ == "__main__":
    main()