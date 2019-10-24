import pandas
import datetime
import os
import shutil
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
        files_tuples.append((row["Read1"], row["Read2"], str(row["Samples"])))
    return files_tuples

"""
Trimmomatic call.

Returns:
    int -- Execution state (0 if everything is all right)
"""
def trimmomatic_call(input_file1, input_file2, phred, trimfile,
                    paired_out_file1, paired_out_file2, unpaired_out_file1, unpaired_out_file2):
    arguments = ["trimmomatic", "PE", phred, input_file1, input_file2, \
                paired_out_file1, paired_out_file2, unpaired_out_file1, unpaired_out_file2, trimfile]
    return call(arguments)

"""
Prinseq call

Returns:
    int -- Execution state (0 if everything is all right)
"""
def prinseq_call(input_file1, input_file2, min_len=40, min_qual_mean=25, trim_qual_right=25, trim_qual_window=15, trim_qual_type="mean", out_format=3, log_name=None):
    arguments = ["prinseq-lite.pl", "-verbose", "-fastq", input_file1, "-fastq2", input_file2, "-min_len", min_len, \
                "-min_qual_mean", min_qual_mean, "-trim_qual_right", trim_qual_right, "-trim_qual_window", \
                trim_qual_window, "-trim_qual_type", trim_qual_type, "-out_format", out_format, "-out_bad", "null", "-log", log_name]
    return call(arguments)

"""
Places prinseq output files in directories with the following structure: /OUTPUT[timestamp]/Prinseq_filtering2/sample
"""
def refactor_prinseq_output(input_dir, output_dir):
    for root, dirs, files in os.walk(input_dir):
        main_out_folder = root.split("/")[0]
        for file in files:
            if file.__contains__("prinseq"):
                strain = file.split("_")[0]
                shutil.move(os.path.join(root, file), main_out_folder+"/"+output_dir+"/"+strain)

def spades_call():
    pass

def contigs_rename_short():
    pass

def quast_call():
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
    output_folder = "Workflow_OUTPUT_"+str(now.strftime('%Y%m%d_%H%M%S'))
    os.mkdir(output_folder)
    trimmomatic_dir = "Trimmomatic_filtering1"
    prinseq_dir = "Prinseq_filtering2"
    os.mkdir(output_folder+"/"+trimmomatic_dir)
    os.mkdir(output_folder+"/"+prinseq_dir)

    # Calls trimmomatic for every row in input_files.csv
    for sample1, sample2, sample_basename in read_input_files("input_files.csv"):
        trimmomatic_call(input_file1=sample1,
                        input_file2=sample2,
                        phred="-phred33",
                        trimfile="ILLUMINACLIP:adapters.fa:1:30:11",
                        paired_out_file1=output_folder+"/"+trimmomatic_dir+"/"+sample_basename+"_R1_paired.fastq",
                        paired_out_file2=output_folder+"/"+trimmomatic_dir+"/"+sample_basename+"_R2_paired.fastq",
                        unpaired_out_file1=output_folder+"/"+trimmomatic_dir+"/"+sample_basename+"_R1_unpaired.fastq",
                        unpaired_out_file2=output_folder+"/"+trimmomatic_dir+"/"+sample_basename+"_R2_unpaired.fastq")
        
        os.mkdir(output_folder+"/"+prinseq_dir+"/"+sample_basename)
        prinseq_call(input_file1=output_folder+"/"+trimmomatic_dir+"/"+sample_basename+"_R1_paired.fastq",
                     input_file2=output_folder+"/"+trimmomatic_dir+"/"+sample_basename+"_R2_paired.fastq", 
                     min_len="40", 
                     min_qual_mean="25", 
                     trim_qual_right="25", 
                     trim_qual_window="15", 
                     trim_qual_type="mean",
                     out_format="3",
                     log_name=output_folder+"/"+prinseq_dir+"/"+sample_basename+"/"+sample_basename+".log")

    refactor_prinseq_output(output_folder+"/"+trimmomatic_dir, prinseq_dir)
    spades_call()
    contigs_rename_short()
    mlst_call()
    abricate_call()
    prokka_call()
    roary_call()