import os
import sqlite3
from Bio import SeqIO
from Levenshtein import distance
import subprocess
import pandas as pd
import numpy as np
import statistics
import shutil

from task1 import *
from utils import *



def filter_sequences(conn,cursor, fastq_name,input_fastq1, input_fastq2, threshold, outdir,min_overlap, min_amplicon_length,fixed_sequence1=None, fixed_sequence2=None):
    """
    Filter sequences from paired-end FASTQ files and update the database.

    Args:
        fastq_name (str): Name of the FASTQ file.
        input_fastq1 (str): Path to the first input FASTQ file.
        input_fastq2 (str): Path to the second input FASTQ file.
        threshold (int): Threshold value for sequence similarity.
        output_fasta (str): Path to the output FASTA file.
        outdir (str): Directory to store output files.
        min_overlap(int): Minimum overlap length 
        min_amplicon_length(int): Minimum Amplicon length
        fixed_sequence1 (str, optional): Fixed sequence for input FASTQ 1. Defaults to None.
        fixed_sequence2 (str, optional): Fixed sequence for input FASTQ 2. Defaults to None.
    """
   
    # 1. merge fastqs in tmp folders
    newoutdir= outdir +"/tmp"
    if not os.path.exists(newoutdir):
        # Create the directory
        os.mkdir(newoutdir)
    processed_fastq= processed_paired_end_fastq(input_fastq1, input_fastq2, min_overlap, min_amplicon_length,newoutdir, fixed_sequence1, fixed_sequence2)
    
  
   
    cursor.execute(f"SELECT * FROM fastq_files WHERE fastq_name = ?", (fastq_name,)) ## wird'es in db vll runtime problematisch
    results= cursor.fetchall()
    
    if len(results)> 0: 
        print("the fastq file is already found in the database. Stopping...")
        return
    

    df = pd.DataFrame(columns=['seq_id', 'intra_freq', 'avg_inter_freq'])
    
    unique_sequences_df = pd.DataFrame(columns=['sequence', 'number_seq'])
    for record in SeqIO.parse(processed_fastq, "fastq"): # multiprocessing (?)
        sequence = str(record.seq)
        #print(sequence)
        df, unique_sequences_df=checkEntry(cursor, sequence, df, unique_sequences_df,threshold)

        print(df)
        print(unique_sequences_df)
                
                

    # create fasta
    # Write unique sequences to FASTA file
    output_fasta=outdir+"/output.fasta"
    with open(output_fasta, 'w') as f:
        for idx, seq in enumerate(unique_sequences_df['sequence'], start=1):
            f.write(f'>Sequence_{idx}\n{seq}\n')
    


     # create csv
    df.to_csv(outdir+'/output.csv', index=False)


    # Update database with new sequence
    updateNewEntry(conn, cursor, fastq_name, df, unique_sequences_df)
    # delete tmp folder
    if os.path.exists(newoutdir):
        shutil.rmtree(newoutdir)
        print(f"Folder '{newoutdir}' deleted successfully.")
    
    conn.close()
    

