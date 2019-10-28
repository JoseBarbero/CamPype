import pandas
import datetime
import os
import shutil
from subprocess import call
from Bio import SeqIO


def read_input_files(indexfile):
    """
    Gets every pair of reads on input_files.csv

    Returns:
        list of tuples -- List of tuples each containing forward read file path, reverse read file path and file basename (just the sample number)
    """
    files_df = pandas.read_csv(indexfile, sep="\t")
    files_tuples = []
    for _, row in files_df.iterrows():
        files_tuples.append((row["Read1"], row["Read2"], str(row["Samples"])))

    return files_tuples


def trimmomatic_call(input_file1, input_file2, phred, trimfile,
                    paired_out_file1, paired_out_file2, unpaired_out_file1, unpaired_out_file2):
    """
    Trimmomatic call.

    Returns:
        int -- Execution state (0 if everything is all right)
    """
    arguments = ["trimmomatic", "PE", phred, input_file1, input_file2, \
                paired_out_file1, unpaired_out_file1, paired_out_file2, unpaired_out_file2, trimfile]
    return call(arguments)


def prinseq_call(input_file1, input_file2, min_len=40, min_qual_mean=25, trim_qual_right=25, trim_qual_window=15, trim_qual_type="mean", out_format=3, log_name=None):
    """
    Prinseq call

    Returns:
        int -- Execution state (0 if everything is all right)
    """
    arguments = ["prinseq-lite.pl", "-verbose", "-fastq", input_file1, "-fastq2", input_file2, "-min_len", min_len, \
                "-min_qual_mean", min_qual_mean, "-trim_qual_right", trim_qual_right, "-trim_qual_window", \
                trim_qual_window, "-trim_qual_type", trim_qual_type, "-out_format", out_format, "-out_bad", "null", "-log", log_name]
    return call(arguments)


def refactor_prinseq_output(input_dir, output_dir, sample):
    """
    Places prinseq output files into directories with the following structure: /OUTPUT[timestamp]/Prinseq_filtering2/sample

    Returns:
        dict -- names of refactored files. key: forward or reverse (R1 or R2), value: filename
    """
    filenames = dict()  # Files with good sequences (except singletons)
    for root, dirs, files in os.walk(input_dir):
        main_out_folder = root.split("/")[0]
        for filename in files:
            if filename.__contains__("prinseq"):
                shutil.move(os.path.join(root, filename), main_out_folder+"/"+output_dir+"/"+sample)
                if filename.startswith(sample+"_R1") and not filename.__contains__("singletons"):
                    filenames["R1"] = filename
                elif filename.startswith(sample+"_R2") and not filename.__contains__("singletons"):
                    filenames["R2"] = filename
    return filenames


def spades_call(forward_sample, reverse_sample, sample, out_dir):
    """
    Spades call

    Returns:
        int -- Execution state (0 if everything is all right)
    """
    arguments = ["spades.py", "-1", forward_sample, "-2", reverse_sample, "--careful", "-o", out_dir+"/"+sample]
    return call(arguments)


def contigs_trim_and_rename(contigs_file, output_dir, min_len):
    """
    Creates new fasta file filtering sequences shorter than min_len and shortening sequence identifiers
    """
    large_sequences = []
    for record in SeqIO.parse(contigs_file, "fasta"):
        if len(record.seq) > min_len:
            record.id = "C_"+"_".join(record.id.split("_")[1:4])
            record.description = ""
            large_sequences.append(record)
    SeqIO.write(large_sequences, output_dir, "fasta")


def quast_call(input_file, output_file, min_contig):
    arguments = ["quast", input_file, "-o", output_file, "--min-contig", str(min_contig), "--no-icarus", "--silent"]
    return call(arguments)

def mlst_call():
    pass

def abricate_call():
    pass

def prokka_call():
    pass

def roary_call():
    pass

if __name__ == "__main__":

    # Create output directories
    now = datetime.datetime.now()
    
    output_folder = "Workflow_OUTPUT_"+str(now.strftime('%Y%m%d_%H%M%S'))
    trimmomatic_dir = "Trimmomatic_filtering1"
    prinseq_dir = "Prinseq_filtering2"
    spades_dir = "SPAdes_assembly"
    contigs_dir = "Contigs_renamed_shorten"
    quast_dir = "Sample_assembly_statistics"

    os.mkdir(output_folder)
    os.mkdir(output_folder+"/"+trimmomatic_dir)
    os.mkdir(output_folder+"/"+prinseq_dir)
    os.mkdir(output_folder+"/"+spades_dir)
    os.mkdir(output_folder+"/"+contigs_dir)
    

    for sample1, sample2, sample_basename in read_input_files("input_files.csv"):
        # Trimmomatic call
        trimmomatic_call(input_file1=sample1,
                        input_file2=sample2,
                        phred="-phred33",
                        trimfile="ILLUMINACLIP:adapters.fa:1:30:11",
                        paired_out_file1=output_folder+"/"+trimmomatic_dir+"/"+sample_basename+"_R1_paired.fastq",
                        unpaired_out_file1=output_folder+"/"+trimmomatic_dir+"/"+sample_basename+"_R1_unpaired.fastq",
                        paired_out_file2=output_folder+"/"+trimmomatic_dir+"/"+sample_basename+"_R2_paired.fastq",
                        unpaired_out_file2=output_folder+"/"+trimmomatic_dir+"/"+sample_basename+"_R2_unpaired.fastq")

        # Creates prinseq output directories
        os.mkdir(output_folder+"/"+prinseq_dir+"/"+sample_basename)

        # Prinseq call
        prinseq_call(input_file1=output_folder+"/"+trimmomatic_dir+"/"+sample_basename+"_R1_paired.fastq",
                     input_file2=output_folder+"/"+trimmomatic_dir+"/"+sample_basename+"_R2_paired.fastq", 
                     min_len="40", 
                     min_qual_mean="25", 
                     trim_qual_right="25", 
                     trim_qual_window="15", 
                     trim_qual_type="mean",
                     out_format="3",
                     log_name=output_folder+"/"+prinseq_dir+"/"+sample_basename+"/"+sample_basename+".log")

        # Prinseq output files refactor
        prinseq_files = refactor_prinseq_output(output_folder+"/"+trimmomatic_dir, prinseq_dir, sample_basename)
        
        # Creates SPAdes output directories
        os.mkdir(output_folder+"/"+spades_dir+"/"+sample_basename)

        # SPAdes call
        spades_call(forward_sample=output_folder+"/"+prinseq_dir+"/"+sample_basename+"/"+prinseq_files["R1"],
                    reverse_sample=output_folder+"/"+prinseq_dir+"/"+sample_basename+"/"+prinseq_files["R2"],
                    sample=sample_basename,
                    out_dir=output_folder+"/"+spades_dir)

        # Trim short contigs and shorten sequences id
        contigs_trim_and_rename(output_folder+"/"+spades_dir+"/"+sample_basename+"/"+"contigs.fasta", 
                                output_folder+"/"+contigs_dir+"/"+sample_basename+"_contigs.fasta",
                                200)

        # Creates Quast output directories

        os.mkdir(output_folder+"/"+spades_dir+"/"+sample_basename+"/"+quast_dir)

        # Quast call
        quast_call( input_file=output_folder+"/"+contigs_dir+"/"+sample_basename+"_contigs.fasta",
                    output_file=output_folder+"/"+spades_dir+"/"+sample_basename+"/"+quast_dir,
                    min_contig=200)


    mlst_call()
    abricate_call()
    prokka_call()
    roary_call()