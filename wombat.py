import pandas
import datetime
import os
from subprocess import call

"""
Gets every pair of reads on input_files.csv

Returns:
    list of tuples -- List of tuples each containing forward read file path, reverse read file path and file basename (just the sample number)
"""
def read_input_files(indexfile):
    files_df = pandas.read_csv(indexfile, sep="\t")
    files_tuples = []
    for _, row in files_df.iterrows():
        sample_basename = row["Read1"].split("/")[-1].split("_")[0]
        files_tuples.append((row["Read1"], row["Read2"], sample_basename))
    return files_tuples

def trimmomatic_call(inputFile1, inputFile2, phred="-phred33", pairedOutputFile1=None, pairedOutputFile2=None, unpairedOutputFile1=None, unpairedOutputFile2=None, trimfile=None):
    arguments = ["trimmomatic", "PE", phred, inputFile1, inputFile2, pairedOutputFile1, pairedOutputFile2, unpairedOutputFile1, unpairedOutputFile2, trimfile]
    
    call(arguments)

def prinseq_call():
    pass

def spades_call():
    pass

def contigs_renamed_shorten():
    pass

def mlst_call():
    pass

def abricate_call():
    pass

def prokka_call():
    pass

def roary_call():
    pass

if __name__ == "__main__":

    # Creates output directory with a timestamp
    now = datetime.datetime.now()
    output_folder = "Workflow_OUTPUT_"+str(now.strftime('%Y%m%d_%H%M')+"/")
    os.mkdir(output_folder)

    # Calls trimmomatic for every row in input_files.csv
    for sample1, sample2, sample_basename in read_input_files("input_files.csv"):
        trimmomatic_call(inputFile1=sample1,
                        inputFile2=sample2,
                        pairedOutputFile1=output_folder+sample_basename+"_R1_paired.fastq",
                        pairedOutputFile2=output_folder+sample_basename+"_R2_paired.fastq",
                        unpairedOutputFile1=output_folder+sample_basename+"_R1_unpaired.fastq",
                        unpairedOutputFile2=output_folder+sample_basename+"_R2_unpaired.fastq",
                        trimfile="ILLUMINACLIP:adapters.fa:1:30:11")

    prinseq_call()
    spades_call()
    contigs_renamed_shorten()
    mlst_call()
    abricate_call()
    prokka_call()
    roary_call()