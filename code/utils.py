import zlib
import glob
import os
import subprocess


dir_ngmerge="./Tools/NGmerge"


def processed_paired_end_fastq(forward_fastq, reverse_fastq,  min_overlap, min_amplicon_length,outdir,fixed_sequence1=None, fixed_sequence2=None):
    """
    Process paired-end FASTQ files by removing adapters and merging overlapping reads.

    Args:
        forward_fastq (str): Path to the forward FASTQ file.
        reverse_fastq (str): Path to the reverse FASTQ file.
        min_overlap (int): Minimum overlap length for merging reads.
        min_amplicon_length (int): Minimum length of amplicon after filtering.
        outdir (str): Directory to save output files.
        fixed_sequence1 (str, optional): Adapter sequence for the forward read. Defaults to None.
        fixed_sequence2 (str, optional): Adapter sequence for the reverse read. Defaults to None.

    Returns:
        str: Path to the processed FASTQ file.
    """
    output_forward_paired=outdir+"/forward.fastq"
    output_reverse_paired=outdir+"/reverse.fastq"
    merge_paired=outdir+"/merge.fastq"
    end_res=outdir+"/end_res.fastq"
    
    ##remove adapter (fixed_sequences)
    if fixed_sequence1 and fixed_sequence2:
        trimm_command= [
        "cutadapt",
        "-a", fixed_sequence1, 
        "-A", fixed_sequence2,
        "-o", output_forward_paired,
        "-p", output_reverse_paired,
        forward_fastq, reverse_fastq
        ]
    else: 
        trimm_command= [
        "cutadapt",
        "-o", output_forward_paired,
        "-p", output_reverse_paired,
        forward_fastq, reverse_fastq
        ]
    subprocess.run(trimm_command)


    ##merging
    merge_command = [
        dir_ngmerge, 
        '-1', output_forward_paired, 
        '-2', output_reverse_paired, 
        '-o', merge_paired,
        '-m', str(min_overlap),
        '-u', '41',
        '-g'
    ]


    subprocess.run(merge_command)

    

    ##filtering by min length
    len_command = ["seqtk", "seq", "-L", str(min_amplicon_length), merge_paired]

    # Execute the command and redirect output to the output file
    with open(end_res, "w") as out_handle:
        subprocess.run(len_command, stdout=out_handle, check=True)

    return end_res


def get_pair_end_fastq_files(directory):
    """
    Get unique pair-end FASTQ file names from a directory.

    Args:
        directory (str): Path to the directory containing FASTQ files.

    Returns:
        set: Set of unique pair-end FASTQ file names.
    """
    fastq_files = glob.glob(directory + "/*.fastq.gz")
    pair_end_files = [os.path.basename(file).split("_")[0] for file in fastq_files if "_1.fastq.gz" in file or "_2.fastq.gz" in file]
    pair_end_files = set(pair_end_files)
    return pair_end_files

def compress_text(text):
    """
    Compress text using zlib compression.

    Args:
        text (str): Text to be compressed.

    Returns:
        bytes: Compressed text.
    """
    # Convert text to bytes and compress it using zlib
    compressed_text = zlib.compress(text.encode())
    return compressed_text

def decompress_text(compressed_text):
    """
    Decompress text using zlib decompression.

    Args:
        compressed_text (bytes): Compressed text.

    Returns:
        str: Decompressed text.
    """
    # Decompress the compressed text using zlib and convert it back to text
    decompressed_text = zlib.decompress(compressed_text).decode()
    return decompressed_text
