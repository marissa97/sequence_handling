import sqlite3
from Bio import SeqIO
from Levenshtein import distance
import pandas as pd
from utils import *
import statistics
import shutil

def activateDB(db_name):
    """
    Connect to the SQLite database and activate necessary functions.

    Args:
        db_name (str): Name of the SQLite database.

    Returns:
        tuple: A tuple containing the connection object and cursor object.
    """
    conn = sqlite3.connect(db_name)
    addEditCalcSQL(conn)
    cursor = conn.cursor()
    return conn, cursor

def addEditCalcSQL(conn):
    """
    Add custom SQL functions to the SQLite connection.

    Args:
        conn: SQLite connection object.
    """
    conn.create_function("levenshtein", 2, distance)
    conn.create_function("decompress_text", 1, decompress_text)

def createTables(conn,cursor):
    """
    Create necessary tables in the SQLite database if they don't exist.

    Args:
        conn: SQLite connection object.
        cursor: SQLite cursor object.
    """
    # Create the first table: sequences
    cursor.execute('''CREATE TABLE IF NOT EXISTS sequences (
                        seq_id INTEGER PRIMARY KEY,
                        sequence TEXT NOT NULL
                    )''')

    # Create the second table: fastq_files
    cursor.execute('''CREATE TABLE IF NOT EXISTS fastq_files (
                        fastq_id INTEGER PRIMARY KEY,
                        fastq_name TEXT
                    )''')
    
    # Create the third table: fastq_sequences
    cursor.execute('''CREATE TABLE IF NOT EXISTS fastq_sequences (
                    rel_id INTEGER PRIMARY KEY,
                    fastq_id INTEGER NOT NULL,
                    seq_id INTEGER NOT NULL,
                    number_seq INTEGER,
                    FOREIGN KEY (seq_id) REFERENCES sequences(seq_id),
                    FOREIGN KEY (fastq_id) REFERENCES fastq_files(fastq_id)
                )''')

    

    # Commit changes and close the connection
    conn.commit()
    
    print("Database created successfully.")



def insert_toDB(conn, cursor, table_name, values):
    """
    Insert values into the specified table in the database.

    Args:
        conn: SQLite connection object.
        cursor: SQLite cursor object.
        table_name (str): Name of the table.
        values: Values to be inserted into the table.
    """

    # Create the INSERT statement dynamically
    placeholders = ', '.join('?' * len(values))
    query = f"INSERT OR IGNORE INTO {table_name} VALUES (NULL,{placeholders})"

    # Execute the INSERT statement with the given values
    cursor.execute(query, values)
    conn.commit()


def select_fromDB(cursor, table_name, toselect, condition_columns=None, condition_symbols=None, condition_values=None):
    """
    Select data from the specified table in the database.

    Args:
        cursor: SQLite cursor object.
        table_name (str): Name of the table.
        toselect (str): Columns to select.
        condition_columns (str): Columns to filter on.
        condition_symbols (str): Symbols for the condition.
        condition_values: Values for the condition.

    Returns:
        list: List of selected rows.
    """

    # Construct the SELECT query dynamically
    query = f"SELECT {toselect} FROM {table_name}"
    if condition_columns and condition_values:
        print(condition_values)
        placeholders = ', '.join('?' * len(condition_values))
        query += f" WHERE {condition_columns} {condition_symbols} {placeholders}"
        print(query)

    
    # Execute the SELECT query with the given condition
    if condition_columns and condition_values:
        cursor.execute(query, (condition_values,))
    else:
        cursor.execute(query)

    # Fetch the results
    rows = cursor.fetchall()

    return rows


def checkExistFastqFile(cursor, fastq_name):
    """
    Check if a FASTQ file exists in the database.

    Args:
        cursor: SQLite cursor object.
        fastq_name (str): Name of the FASTQ file.

    Returns:
        list: List of matching rows.
    """
    cursor.execute(f"SELECT * FROM fastq_files WHERE fastq_name = ?", (fastq_name,)) ## wird'es in db vll runtime problematisch
    results= cursor.fetchall()
    return results

def checkEntry(cursor, sequence, df, unique_sequences_df, threshold):
    """
    Check if a sequence exists in the database, update dataframes accordingly.

    Args:
        cursor: SQLite cursor object.
        sequence (str): Sequence to check.
        df: Dataframe containing sequence information.
        unique_sequences_df: Dataframe containing unique sequences and their counts.
        threshold (int): Threshold value for Levenshtein distance.

    Returns:
        tuple: Updated dataframes.
    """
    
    sim_id=None

    ## select exacte seq in db 
    cursor.execute(f"SELECT seq_id FROM sequences WHERE decompress_text(sequence) = ?", (sequence,)) ## wird'es in db vll runtime problematisch
    results= cursor.fetchone()
    #print(results)
    sim_id=None
    # if false  -> levenshtein, select one (runtime)
    if results is not None:
        sim_id= results[0]
                
                
    else: 
        
        cursor.execute(f"SELECT seq_id FROM sequences WHERE levenshtein(decompress_text(sequence),  ?) < ? LIMIT 1", (sequence, threshold,)) ## wird'es in db vll runtime problematisch
        results = cursor.fetchone()
        print("levenshtein:", results)
        if results is not None: 
            sim_id=results[0]

    if sim_id is not None:
        if sim_id not in df['seq_id']: ## similarity in df??
            query= f"SELECT * FROM fastq_sequences WHERE seq_id = {sim_id}" # nur 1 seq
            cursor.execute(query)
            results = cursor.fetchall()
            avg_inter_freq= statistics.mean([row[3] for row in results])
            selected_df= pd.DataFrame({'seq_id': [sim_id], 'avg_inter_freq': [avg_inter_freq]})
                

            
                    
                
            # if sim_id already exist in df?
            # select_seqID, already in df
            # Merge dataframes on seq_id with an outer join
            merged_df = pd.merge(df, selected_df, on='seq_id', how='outer', suffixes=('_df', '_selected'))

            # Filter out rows where seq_id is present in both dataframes
            new_rows = merged_df[merged_df['avg_inter_freq_selected'].notna() & merged_df['avg_inter_freq_df'].isna()]

                    # Concatenate the remaining rows with df
            df = pd.concat([df, new_rows[['seq_id', 'avg_inter_freq_selected']].rename(columns={'avg_inter_freq_selected': 'avg_inter_freq'})], ignore_index=True)

            df.loc[df['seq_id']== sim_id, 'intra_freq']=1
                
        else:
            # add intra_frequency +1    
            df.loc[df['seq_id']== sim_id, 'intra_freq'] = df.loc[df['seq_id']== sim_id, 'intra_freq']+1 

    #if not similar_sequence:
    else: 
        if sequence not in unique_sequences_df['sequence'].values:
            unique_sequences_df = unique_sequences_df._append({'sequence': sequence, 'number_seq': 1}, ignore_index=True)
        else: 
            unique_sequences_df.loc[unique_sequences_df['sequence'] == sequence, 'number_seq'] = unique_sequences_df.loc[unique_sequences_df['sequence'] == sequence, 'number_seq']+1 

    return df, unique_sequences_df               
                
              

         
            
def updateNewEntry(conn, cursor, fastq_name, df, unique_sequences_df):
    """
    Update database with new sequences.

    Args:
        conn: SQLite connection object.
        cursor: SQLite cursor object.
        fastq_name (str): Name of the FASTQ file.
        df: Dataframe containing sequence information.
        unique_sequences_df: Dataframe containing unique sequences and their counts.
    """
    print("updating db")
    # Update database with new sequence
    # insert into 'sequences'
    unique_sequences_df['seq_id'] = None
    seq_df= unique_sequences_df[['seq_id', 'sequence']]
    seq_df['sequence']= [compress_text(x) for x in seq_df['sequence']]
    seq_df.to_sql('sequences', conn, if_exists='append', index=False)
    #get seq_id
    cursor.execute("SELECT last_insert_rowid()")
    
    seq_id_last = cursor.fetchone()[0]
    print(unique_sequences_df)
    print("seq_id: ",seq_id_last)

    
    # insert into 'fastq_files'
    query = f"INSERT OR IGNORE INTO fastq_files VALUES (NULL,?)"

    # Execute the INSERT statement with the given values
    cursor.execute(query, (fastq_name,))
    conn.commit()
    fastq_id_last = cursor.lastrowid
    print(fastq_id_last)

    # insert into fastq_sequences (rel_id, seq_id, fastq_id, number_seq)
    
    seq_df['seq_id']= range(seq_id_last-len(seq_df)+1, seq_id_last+1)
    seq_df['fastq_id']= fastq_id_last
    seq_df['rel_id']= None
    fastq_seq_df= seq_df[['rel_id', 'fastq_id', 'seq_id']]
    fastq_seq_df['number_seq']=unique_sequences_df['number_seq']
    fastq_seq_df.to_sql('fastq_sequences', conn, if_exists='append', index=False)

    # Update database with found sequences in filtered fastq files
    # sim_id, fastq_id, seq_id:None, number_seq
    found_df= pd.DataFrame({'rel_id': None, 'fastq_id': fastq_id_last, 'seq_id': df['seq_id'], 'number_seq': df['intra_freq']})
    found_df.to_sql('fastq_sequences', conn, if_exists='append', index=False)
    
    

def mainCreateDB(conn, cursor, inputFolder, min_overlap, min_amplicon_length,outdir, threshold, fixed_sequence1=None, fixed_sequence2=None):
    """
    Main function to create and update the database with information from FASTQ files.

    Args:
        conn: SQLite connection object.
        cursor: SQLite cursor object.
        inputFolder (str): Path to the folder containing FASTQ files.
        min_overlap (int): Minimum overlap required for sequence processing.
        min_amplicon_length (int): Minimum amplicon length for sequence processing.
        outdir (str): Directory to store processed files.
        threshold (int): Threshold value for sequence similarity.
        fixed_sequence1 (str, optional): Fixed sequence for input FASTQ 1. Defaults to None.
        fixed_sequence2 (str, optional): Fixed sequence for input FASTQ 2. Defaults to None.
    """

    createTables(conn, cursor)
    print(inputFolder)
    fastq_files_name= get_pair_end_fastq_files(inputFolder)
    

    for fastq_file in fastq_files_name:
        results= checkExistFastqFile(cursor, fastq_file)
        
        if len(results) >0: 
            print("the fastq file is already found in the database.")
        
        else:
            #get input1 and input2:
            print(inputFolder)
            input1= inputFolder+"/"+fastq_file+"_1.fastq.gz"
            input2= inputFolder+"/"+fastq_file+"_2.fastq.gz"


            df = pd.DataFrame(columns=['seq_id', 'intra_freq', 'avg_inter_freq'])
            unique_sequences_df = pd.DataFrame(columns=['sequence', 'number_seq'])

            #processing fastq_files-> tmp directory for fastq files processing
            newoutdir= outdir +"/"+ fastq_file
            if not os.path.exists(newoutdir):
                # Create the directory
                os.mkdir(newoutdir)
            processed_fastq=processed_paired_end_fastq(input1, input2, min_overlap, min_amplicon_length,newoutdir, fixed_sequence1, fixed_sequence2)
            
            #check and update
            for record in SeqIO.parse(processed_fastq, "fastq"): # multiprocessing (?)
                sequence = str(record.seq)
                #print(sequence)
                df, unique_sequences_df=checkEntry(cursor, sequence, df, unique_sequences_df,threshold)

            print(df)
            print(unique_sequences_df)
        
            updateNewEntry(conn, cursor, fastq_file, df, unique_sequences_df)
            print("updated")

            # delete tmp folder
            if os.path.exists(newoutdir):
                shutil.rmtree(newoutdir)
                print(f"Folder '{newoutdir}' deleted successfully.")
    





