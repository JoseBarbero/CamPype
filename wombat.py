import pandas
import datetime
import os
import shutil
import logging
import sys
from terminal_banner import Banner
from subprocess import call
from Bio import SeqIO


def read_aux_files(inputfile):
    """
    File containing auxiliary files routes.
    
    Arguments:
        inputfile {string} -- Filename (and route) of the file containing inputs.
    
    Returns:
        files_routes {list} -- List of auxiliary files.
    """
    files_df = pandas.read_csv(inputfile, sep="\t")
    files_routes = []
    for _, row in files_df.iterrows():
        files_routes.append(row["Route"])
    return files_routes


def read_input_files(indexfile):
    """
    Gets every pair of reads on input_files.csv
    
    Arguments:
        indexfile {string} -- Filename (and route) of the file containing inputs.
    
    Returns:
        files_tuples {list of tuples} -- List of tuples each containing forward read file path, reverse read file path and file basename (just the sample number).
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
    
    Arguments:
        input_file1 {string} -- Input file forward (and route).
        input_file2 {string} -- Input file reverse (and route).
        phred {string} -- Trimmomatic phred parameter.
        trimfile {string} -- File with trimming sequences.
        paired_out_file1 {string} -- Forward file paired output.
        paired_out_file2 {string} -- Reverse file paired output.
        unpaired_out_file1 {string} -- Forward file unpaired output.
        unpaired_out_file2 {string} -- Reverse file unpaired output.
    
    Returns:
        {int} -- Execution state (0 if everything is all right)
    """
    arguments = ["trimmomatic", "PE", phred, input_file1, input_file2, \
                paired_out_file1, unpaired_out_file1, paired_out_file2, unpaired_out_file2, trimfile]
    return call(arguments)


def prinseq_call(input_file1, input_file2, min_len=40, min_qual_mean=25, trim_qual_right=25,
                 trim_qual_window=15, trim_qual_type="mean", out_format=3, log_name=None):
    """
    Prinseq call
    
    Arguments:
        input_file1 {string} -- Input file forward (and route).
        input_file2 {string} -- Input file reverse (and route).
    
    Keyword Arguments:
        min_len {int} -- Minimum read length. (default: {40})
        min_qual_mean {int} -- Minimum read quality. (default: {25})
        trim_qual_right {int} -- Trim sequence by quality score from the 3'-end with this threshold score. (default: {25})
        trim_qual_window {int} -- Trim sequence by quality score from the 5'-end with this threshold score. (default: {15})
        trim_qual_type {str} -- Type of quality score calculation to use. (default: {"mean"})
        out_format {int} -- Output format 1 (FASTA only), 2 (FASTA and QUAL), 3 (FASTQ), 4 (FASTQ and FASTA), 5 (FASTQ, FASTA and QUAL) (default: {3})
        log_name {string} -- Output log file name.
    
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
    
    Arguments:
        input_dir {string} -- Input directory.
        output_dir {string} -- Output directory.
        sample {string} -- Sample basename.
    
    Returns:
        {dict} -- names of refactored files. key: forward or reverse (R1 or R2), value: filename
    """
    filenames = dict()  # Files with good sequences (except singletons)
    for root, _dirs, files in os.walk(input_dir):
        main_out_folder = root.split("/")[0]
        for filename in files:
            if filename.__contains__("prinseq"):
                if filename.__contains__("singletons"):
                    os.remove(os.path.join(root, filename))
                else:
                    # Move every prinseq file from trimmomatic folder to prinseq folder
                    shutil.move(os.path.join(root, filename), main_out_folder+"/"+output_dir+"/"+sample)
                    if filename.startswith(sample+"_R1"):
                        filenames["R1"] = filename
                    elif filename.startswith(sample+"_R2"):
                        filenames["R2"] = filename
    return filenames


def spades_call(forward_sample, reverse_sample, sample, out_dir):
    """
    Spades call
    
    Arguments:
        forward_sample {string} -- Forward sample file name (and route).
        reverse_sample {string} -- Reverse sample file name (and route).
        sample {string} -- Sample basename.
        out_dir {string} -- Output directory.
    
    Returns:
        {int} -- Execution state (0 if everything is all right)
    """
    arguments = ["spades.py", "-1", forward_sample, "-2", reverse_sample, "--careful", "--cov-cutoff", "auto", "-o", out_dir+"/"+sample]    
    return call(arguments)


def contigs_trim_and_rename(contigs_file, output_dir, min_len):
    """
    Creates new fasta file filtering sequences shorter than min_len and shortening sequence identifiers.
    
    Arguments:
        contigs_file {string} -- Original contigs filename (and route).
        output_dir {string} -- Output directory.
        min_len {int} -- Minimum sequence length.
    """
    large_sequences = []
    for record in SeqIO.parse(contigs_file, "fasta"):
        if len(record.seq) > min_len:
            record.id = "C_"+"_".join(record.id.split("_")[1:4])
            record.description = ""
            large_sequences.append(record)
    SeqIO.write(large_sequences, output_dir, "fasta")


def quast_call(input_file, output_dir, min_contig):
    """
    Quast call.
    
    Arguments:
        input_file {string} -- Input file (and route).
        output_dir {string} -- Output directory.
        min_contig {int} -- Lower threshold for a contig length (in bp).
    
    Returns:
        {int} -- Execution state (0 if everything is all right)
    """
    arguments = ["quast", input_file, "-o", output_dir, "--min-contig", str(min_contig), "--no-icarus", "--silent"]
    return call(arguments)


def mlst_call(input_dir, output_dir, output_filename):
    """
    MLST call for every fasta file in input_dir.
    
    Arguments:
        input_dir {string} -- Input directory containing contig files.
        output_dir {string} -- Output directory.
        output_filename {string} -- Output file name (and route).
    
    Returns:
        {int} -- Execution state (0 if everything is all right)
    """
    output_file = open(output_dir+"/"+output_filename, "w")
    
    input_filenames = []

    for _root, _dirs, files in os.walk(input_dir):
        for filename in files:
            if filename.endswith(".fasta"):
                input_filenames.append(input_dir+"/"+filename)

    arguments = ["mlst", *input_filenames]
    return call(arguments, stdout=output_file)


def abricate_call(input_dir, output_dir, output_filename, database):
    """
    ABRicate call.
    
    Arguments:
        input_dir {string} -- Input directory containing contig files.
        output_dir {string} -- Output directory.
        output_filename {string} -- Output file name (and route).
        database {string} -- Database name.
    
    Returns:
        {int} -- Execution state (0 if everything is all right)
    """
    output_file = open(output_dir+"/"+output_filename, "w")
    
    input_filenames = []

    for _root, _dirs, files in os.walk(input_dir):
        for filename in files:
            if filename.endswith(".fasta"):
                input_filenames.append(input_dir+"/"+filename)

    arguments = ["abricate", *input_filenames, "--db", database]
    return call(arguments, stdout=output_file)


def prokka_call(locus_tag, output_dir, prefix, input_file):
    """
    Prokka call.
    
    Arguments:
        locus_tag {string} -- Locus tag prefix.
        output_dir {string} -- Output directory.
        prefix {string} -- Filename output prefix.
        input_file {string} -- Input filename (and route).
    
    Returns:
        {int} -- Execution state (0 if everything is all right)
    """
    arguments = ["prokka", "--locustag", locus_tag, "--outdir", output_dir, "--prefix", prefix, "--kingdom", "Bacteria", "--gcode", "11", input_file]
    return call(arguments)


def dfast_call(input_file, min_length, out_path, sample_basename):
    """
    Dfast call.
    
    Arguments:
        input_file {string} -- Input filename (and route).
        min_length {int} -- Minimum sequence length.
        out_path {strin} -- Output folder.
    
    Returns:
        {int} -- Execution state (0 if everything is all right)
    """
    arguments = ["dfast", "--genome", input_file, "--minimum_length", str(min_length), "--out", out_path]
    state = call(arguments)

    # Replace default output filenames including string basename
    for root, _dirs, files in os.walk(out_path):
        for filename in files:
            if filename.__contains__("genome"):
                new_filename = filename.replace("genome", sample_basename)
                shutil.move(os.path.join(root, filename), os.path.join(root, new_filename))

    return state


def roary_call(input_files, output_dir):
    """
    Roary call.
    
    Arguments:
        input_files {list} -- GFF files from prokka.
        output_dir {string} -- Output directory.
    
    Returns:
        {int} -- Execution state (0 if everything is all right)
    """
    arguments = ["roary", "-f", output_dir, *input_files]
    ex_state = call(arguments)
    # Set Roary output directory name
    for _root, dirs, _files in os.walk("."):
        for dirname in dirs:
            if dirname.startswith("Roary_pangenome_"):
                os.rename(dirname, "Roary_pangenome")
    return ex_state


def roary_plots_call(input_newick, input_gene_presence_absence, output_dir):
    """
    Roary plots call
    
    Arguments:
        input_newick {string} -- Filename (and route) to newick input file.
        input_gene_presence_absence {[type]} -- Filename (and route) to gene presence/absence input file.
        output_dir {[type]} -- Route to output files.
    
    Returns:
        {int} -- Execution state (0 if everything is all right)
    """
    arguments = ["python", "utils/roary_plots.py", input_newick, input_gene_presence_absence]
    ex_state = call(arguments)
    
    # Roary_plots saves output files in the current directory, so we move them to our own
    for root, _dirs, files in os.walk("."):
        for filename in files:
            if filename.startswith("pangenome_"):
                shutil.move(root+"/"+filename, output_dir+"/"+filename)
    return ex_state


if __name__ == "__main__":

    # Get selected annotation tool (prokka by default)
    if len(sys.argv) > 2:
        if sys.argv[2] == "--annotator":
            annotator = sys.argv[3]
    else:
        annotator = "prokka"

    # Create output directories
    now = datetime.datetime.now()

    adapters_file, reference_annotation_file = read_aux_files("reference_files.csv")
    
    output_folder = sys.argv[1]
    trimmomatic_dir = "Trimmomatic_filtering1"
    prinseq_dir = "Prinseq_filtering2"
    spades_dir = "SPAdes_assembly"
    contigs_dir = "Contigs_renamed_shorten"
    mlst_dir = "MLST"
    abricate_vir_dir = "ABRicate_virulence_genes"
    abricate_abr_dir = "ABRicateAntibioticResistanceGenes"
    prokka_dir = "Prokka_annotation"
    dfast_dir = "Dfast_annotation"
    roary_dir = "Roary_pangenome"
    roary_plots_dir = "Roary_plots"

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
        
    os.mkdir(output_folder+"/"+trimmomatic_dir)
    os.mkdir(output_folder+"/"+prinseq_dir)
    os.mkdir(output_folder+"/"+spades_dir)
    os.mkdir(output_folder+"/"+contigs_dir)
    os.mkdir(output_folder+"/"+mlst_dir)
    os.mkdir(output_folder+"/"+abricate_vir_dir)
    os.mkdir(output_folder+"/"+abricate_abr_dir)
    if annotator == "dfast":
        os.mkdir(output_folder+"/"+dfast_dir)
    else:
        os.mkdir(output_folder+"/"+prokka_dir)

    roary_input_files = []
    for sample1, sample2, sample_basename in read_input_files("input_files.csv"):
        # Trimmomatic call
        print(Banner("\nStep 1 for sequence "+sample_basename+": Trimmomatic\n"), flush=True)
        trimmomatic_call(input_file1=sample1,
                        input_file2=sample2,
                        phred="-phred33",
                        trimfile="ILLUMINACLIP:"+adapters_file+":1:30:11",
                        paired_out_file1=output_folder+"/"+trimmomatic_dir+"/"+sample_basename+"_R1_paired.fastq",
                        unpaired_out_file1=output_folder+"/"+trimmomatic_dir+"/"+sample_basename+"_R1_unpaired.fastq",
                        paired_out_file2=output_folder+"/"+trimmomatic_dir+"/"+sample_basename+"_R2_paired.fastq",
                        unpaired_out_file2=output_folder+"/"+trimmomatic_dir+"/"+sample_basename+"_R2_unpaired.fastq")

        # Create prinseq output directories
        os.mkdir(output_folder+"/"+prinseq_dir+"/"+sample_basename)

        # Prinseq call
        print(Banner("\nStep 2 for sequence "+sample_basename+": Prinseq\n"), flush=True)
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
        
        # Create SPAdes output directories
        os.mkdir(output_folder+"/"+spades_dir+"/"+sample_basename)

        # SPAdes call
        print(Banner("\nStep 3 for sequence "+sample_basename+": SPAdes\n"), flush=True)
        spades_call(forward_sample=output_folder+"/"+prinseq_dir+"/"+sample_basename+"/"+prinseq_files["R1"],
                    reverse_sample=output_folder+"/"+prinseq_dir+"/"+sample_basename+"/"+prinseq_files["R2"],
                    sample=sample_basename,
                    out_dir=output_folder+"/"+spades_dir)

        # Trim short contigs and shorten sequences id
        contigs_trim_and_rename(output_folder+"/"+spades_dir+"/"+sample_basename+"/"+"contigs.fasta", 
                                output_folder+"/"+contigs_dir+"/"+sample_basename+"_contigs.fasta",
                                500)    # TODO Para el futuro, este valor deberá ser el doble de la longitud de una read del archivo .fastq del secuenciador

        # Create Quast output directories
        quast_dir = sample_basename+"_assembly_statistics"
        os.mkdir(output_folder+"/"+spades_dir+"/"+sample_basename+"/"+quast_dir)

        # Quast call
        print(Banner("\nStep 4 for sequence "+sample_basename+": Quast\n"), flush=True)
        quast_call( input_file=output_folder+"/"+contigs_dir+"/"+sample_basename+"_contigs.fasta",
                    output_dir=output_folder+"/"+spades_dir+"/"+sample_basename+"/"+quast_dir,
                    min_contig=200)

        # Annotation (Prokka or dfast)
        if annotator.lower() == "dfast":
            # Dfast call
            annotation_dir = dfast_dir
            print(Banner("\nStep 5 for sequence "+sample_basename+": Dfast\n"), flush=True)
            dfast_call(input_file=output_folder+"/"+contigs_dir+"/"+sample_basename+"_contigs.fasta",
                       min_length=0,
                       out_path=output_folder+"/"+dfast_dir+"/"+sample_basename,
                       sample_basename=sample_basename)
        else:
            # Prokka call
            annotation_dir = prokka_dir
            print(Banner("\nStep 5 for sequence "+sample_basename+": Prokka\n"), flush=True)
            prokka_call(locus_tag=sample_basename+"_L",
                        output_dir=output_folder+"/"+prokka_dir+"/"+sample_basename,
                        prefix=sample_basename,
                        input_file=output_folder+"/"+contigs_dir+"/"+sample_basename+"_contigs.fasta")

        # Set roary input files
        roary_input_files.append(output_folder+"/"+annotation_dir+"/"+sample_basename+"/"+sample_basename+".gff")

    # Annotate reference fasta file 
    if reference_annotation_file:
        reference_annotation_filename = reference_annotation_file.split("/")[-1]
        reference_annotation_basename = reference_annotation_filename.split(".")[-2]
        if annotator.lower() == "dfast":
            # Dfast call
            annotation_dir = dfast_dir
            print(Banner("\nStep 5 for reference sequence: Dfast\n"), flush=True)
            dfast_call(input_file=reference_annotation_file,
                        min_length=0,
                        out_path=output_folder+"/"+dfast_dir+"/"+reference_annotation_basename,
                        sample_basename=reference_annotation_basename)
        else:
            # Prokka call
            annotation_dir = prokka_dir
            print(Banner("\nStep 5 for reference sequence: Prokka\n"), flush=True)
            prokka_call(locus_tag=reference_annotation_basename+"_L",
                        output_dir=output_folder+"/"+prokka_dir+"/"+reference_annotation_basename,
                        prefix=reference_annotation_basename,
                        input_file=reference_annotation_file)

        # Set roary input files
        roary_input_files.append(output_folder+"/"+annotation_dir+"/"+reference_annotation_basename+"/"+reference_annotation_basename+".gff")

    # MLST call
    print(Banner("\nStep 6: MLST\n"), flush=True)
    mlst_call(input_dir=output_folder+"/"+contigs_dir,
            output_dir=output_folder+"/"+mlst_dir,
            output_filename="MLST.txt")

    # ABRicate call (virulence genes)
    print(Banner("\nStep 7: ABRicate (virulence genes)\n"), flush=True)
    abricate_call(input_dir=output_folder+"/"+contigs_dir,
                output_dir=output_folder+"/"+abricate_vir_dir,
                output_filename="SampleVirulenceGenes.tab",
                database = "vfdb")

    # ABRicate call (antibiotic resistance genes)
    print(Banner("\nStep 8: ABRicate (antibiotic resistance genes)\n"), flush=True)
    abricate_call(input_dir=output_folder+"/"+contigs_dir,
                output_dir=output_folder+"/"+abricate_abr_dir,
                output_filename="SampleAntibioticResistanceGenes.tab",
                database = "resfinder")

    # Roary call
    print(Banner("\nStep 9: Roary\n"), flush=True)
    roary_call(input_files=roary_input_files, output_dir=output_folder+"/"+roary_dir)

    # Roary plots call
    os.mkdir(output_folder+"/"+roary_dir+"/"+roary_plots_dir)
    print(Banner("\nStep 10: Roary Plots\n"), flush=True)
    roary_plots_call(input_newick=output_folder+"/"+roary_dir+"/accessory_binary_genes.fa.newick",
                    input_gene_presence_absence=output_folder+"/"+roary_dir+"/gene_presence_absence.csv",
                    output_dir=output_folder+"/"+roary_dir+"/"+roary_plots_dir)

    print(Banner("\nDONE\n"), flush=True)
