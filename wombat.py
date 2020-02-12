import pandas as pd
import datetime
import os
import shutil
import re
import logging
import sys
import workflow_config as cfg
import gffutils
from terminal_banner import Banner
from subprocess import call
from Bio import SeqIO


def read_input_files(indexfile):
    """
    Gets every pair of reads on input_files.csv
    
    Arguments:
        indexfile {string} -- Filename (and route) of the file containing inputs.
    
    Returns:
        files_data {list of tuples} -- List of tuples each containing forward read file path, reverse read file path and file basename (just the sample number).
    """
    files_df = pd.read_csv(indexfile, sep="\t")
    files_data = []
    for _, row in files_df.iterrows():
        files_data.append((str(row["Samples"]), {"FW": row["Forward"], "RV": row["Reverse"], "Genus": row["Genus"], "Species": row["Species"]}))
    return files_data


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


def prinseq_call(input_file1, input_file2, output_folder, sample, log_name=None):
    """
    Prinseq call
    
    Arguments:
        input_file1 {string} -- Input file forward (and route).
        input_file2 {string} -- Input file reverse (and route).
        output_folder {string} -- Output folder.
    
    Keyword Arguments:
        log_name {string} -- Output log file name.
    
    Returns:
        int -- Execution state (0 if everything is all right)
    """
    arguments = ["prinseq-lite.pl", "-verbose", "-fastq", input_file1, "-fastq2", input_file2, "-min_len", str(cfg.config["prinseq"]["min_len"]), \
                "-min_qual_mean", str(cfg.config["prinseq"]["min_qual_mean"]), "-trim_qual_right", str(cfg.config["prinseq"]["trim_qual_right"]), "-trim_qual_window", \
                str(cfg.config["prinseq"]["trim_qual_window"]), "-trim_qual_type", cfg.config["prinseq"]["trim_qual_type"], "-out_format", str(cfg.config["prinseq"]["out_format"]), "-out_good", output_folder+"/"+sample, "-out_bad", cfg.config["prinseq"]["out_bad"], "-log", log_name]
    return call(arguments)


def refactor_prinseq_output(input_dir, sample):
    """
    Rename and remove some files from prinseq.
    
    Arguments:
        input_dir {list} -- Input directory.
        sample {string} -- Sample basename.
    
    Returns:
        {dict} -- names of refactored files. key: forward or reverse (R1 or R2), value: filename
    """
    filenames = dict()  # Files with good sequences (except singletons)
    for root, _dirs, files in os.walk(input_dir):
        for filename in files:
            if filename.__contains__("singletons"):
                os.remove(os.path.join(root, filename))
            else:
                if filename.__contains__("_1"):
                    new_filename = os.path.join(root, sample+"_R1.fastq")
                    os.rename(os.path.join(root, filename), new_filename)
                    filenames["R1"] = new_filename
                elif filename.__contains__("_2"):
                    new_filename = os.path.join(root, sample+"_R2.fastq")
                    os.rename(os.path.join(root, filename), new_filename)
                    filenames["R2"] = new_filename
    return filenames


def flash_call(input_file_1, input_file_2, output_filename, output_dir):
    """Flash call.
    
    Arguments:
        input_file_1 {string} -- Forward prinseq file.
        input_file_2 {string} -- Reverse prinseq file.
        output_filename {string} -- Output filename.
        output_dir {string} -- Output folder.
    
    Returns:
        {int} -- Execution state (0 if everything is all right)
    """
    arguments = ["flash", input_file_1, input_file_2, "-o", output_filename, "-d", output_dir]
    return call(arguments)


def spades_call(merged_sample, forward_sample, reverse_sample, sample, out_dir):
    """
    Spades call
    
    Arguments:
        merged_sample {string} -- Merged sample file name (and route).
        forward_sample {string} -- Forward sample file name (and route).
        reverse_sample {string} -- Reverse sample file name (and route).
        sample {string} -- Sample basename.
        out_dir {string} -- Output directory.
    
    Returns:
        {int} -- Execution state (0 if everything is all right)
    """
    arguments = ["spades.py", "--merged", merged_sample, "-1", forward_sample, "-2", reverse_sample, cfg.config["spades"]["mode"], "--cov-cutoff", cfg.config["spades"]["cov_cutoff"], "-o", out_dir+"/"+sample]    
    return call(arguments)


def mauve_call(output_folder, reference_sequence, input_contigs, sample_basename):
    """Mauve call. (Reordering contigs).

    MauveCM will output a series of folders called alignment1-alignmentX, representing each iteration of the reorder.
    
    Arguments:
        output_folder {string} -- Output folder route.
        reference_sequence {string} -- Reference sequence route.
        input_contigs {string} -- Contigs file route.
        sample_basename -- Sample basename.
    
    Returns:
        {string} -- Mauve reordered contigs file path
    """
    arguments = ["MauveCM", "-output", output_folder, "-ref", reference_sequence, "-draft", input_contigs]    
    call(arguments)
    # Here we take the fasta file from the last iteration folder.
    shutil.copyfile(output_folder+"/"+max(next(os.walk(output_folder))[1])+"/"+sample_basename+".fasta", output_folder+"/../contigs/"+sample_basename+".fasta")
    return output_folder+"/../contigs/"+sample_basename+".fasta"


def snippy_call(reference_genome, contigs, output_dir, prefix):
    """
    Snippy call. (SNP identifier)
    
    Arguments:
        referece_genome {string} -- Reference genome file route.
        contigs {string} -- Contigs file route.
        output_dir {string} -- Output directory.
        prefix {string} -- Sample name as prefix.
    
    Returns:
        {int} -- Execution state (0 if everything is all right)
    """
    arguments = ["snippy", "--ref", reference_genome, "--ctgs", contigs, "--outdir", output_dir, "--prefix", prefix, "--report"]    
    return call(arguments)


def contigs_trim_and_rename(contigs_file, output_filename, output_dir, min_len):
    """
    Creates new fasta file filtering sequences shorter than min_len and shortening sequence identifiers.
    
    Arguments:
        contigs_file {string} -- Original contigs filename (and route).
        output_filename {string} -- Output file name.
        output_dir {string} -- Output directory.
        min_len {int} -- Minimum sequence length.
    """
    with open(output_dir+"/Min_contig_len_"+str(min_len)+".txt", 'w') as f:
        f.write("Minimum contig length was set to "+str(min_len))

    large_sequences = []
    for record in SeqIO.parse(contigs_file, "fasta"):
        if len(record.seq) > min_len:
            record.id = "C_"+"_".join(record.id.split("_")[1:4])
            record.description = ""
            large_sequences.append(record)
    SeqIO.write(large_sequences, output_dir+"/"+output_filename, "fasta")


def get_reads_length(input_file):
    """
    Returns first read length form fastq file.
    
    Arguments:
        input_file {string} -- Fastq file.
    """
    first_record = next(SeqIO.parse(input_file, "fastq"))
    return len(first_record.seq)


def quast_call(input_file, output_dir, min_contig_len):
    """
    Quast call.
    
    Arguments:
        input_file {string} -- Input file (and route).
        output_dir {string} -- Output directory.
        min_contig_len -- Lower threshold for a contig length (in bp).
    
    Returns:
        {int} -- Execution state (0 if everything is all right)
    """
    arguments = ["quast", input_file, "-o", output_dir, "--min-contig", str(min_contig_len), cfg.config["quast"]["icarus"], cfg.config["quast"]["mode"]]
    return call(arguments)


def quast_report_unification(input_dir, samples, output_dir):
    """
    Create a report unifying reports from every sample.
    
    Arguments:
        input_dir {string} -- Input directory.
        samples {list} -- List of names of samples.
        output_dir {string} -- Output directory.
    """
    first_col_df = pd.read_csv(input_dir+"/"+samples[0]+"/"+samples[0]+"_assembly_statistics/report.tsv", sep="\t")
    combined_df = pd.concat([pd.read_csv(input_dir+"/"+sample+"/"+sample+"_assembly_statistics/report.tsv", sep="\t").iloc[:, 1] for sample in samples[1:]], axis=1)
    final_df = pd.concat([first_col_df, combined_df], axis=1)
    final_df.to_csv(output_dir+"/quality_assembly_report.tsv", sep="\t")

 
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


def blast_call(proteins_file_ori, proteins_file_dest, contigs_files_paths, blast_database_output, blast_output_folder, blast_output_name):
    """
    Blast call.
    
    Arguments:
        proteins_file_ori {string} -- Reference proteins file path (origin).
        proteins_file_dest {string} -- Reference proteins file path (destination).
        contigs_files_paths {list} -- List of contig files from SPAdes.
        blast_database_output {string} -- Destination folder to blast database.
        blast_output_folder {string} -- Destination folder to blast output.
        blast_output_name {string} -- Output file name.
    
    Returns:
        {int} -- Execution state (0 if everything is all right)
    """
    # Move VF_custom.txt to ABRicateVirulenceGenes/BLAST_custom_proteins
    shutil.copy(proteins_file_ori, proteins_file_dest)

    # Concat every contig from SPAdes on DNA_database.fna (replacing sequences names)
    with open(blast_database_output, "w") as output_file:
        for contig_file_path in contigs_files_paths:
            for record in SeqIO.parse(contig_file_path, "fasta"):
                strain = contig_file_path.split("/")[-2]
                record.id = record.name = record.description = strain+"_N_"+"_".join(record.id.split("_")[1:])
                SeqIO.write(record, output_file, "fasta")
                
    # Create blast database
    blast_db_path = os.path.dirname(os.path.abspath(blast_database_output))+"/DNA_database"
    call(["makeblastdb", "-in", blast_database_output, "-dbtype", cfg.config["blast"]["dbtype"], 
          "-out", blast_db_path, "-title", "DNA_Database"])

    # Call tblastn
    tblastn_state = call(["tblastn", "-db", blast_db_path, "-query", proteins_file_dest, "-evalue", str(cfg.config["blast"]["evalue"]), "-outfmt", 
                        cfg.config["blast"]["outfmt"],
                        "-out", blast_output_folder+"/"+blast_output_name])

    # Add header to tblastn output
    with open(blast_output_folder+"/"+blast_output_name, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        headers = ["query id", "subject id", "% identity", "alignment length", "mismatches", "gap opens", 
                   "query start", "query end", "subject start", "subject end", "evalue", "score", "aligned part of subject sequence"]

        f.write("\t".join(headers) + '\n' + content)

    return tblastn_state


def prokka_call(locus_tag, output_dir, prefix, input_file, genus, species, strain):
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
    arguments = ["prokka", 
                "--locustag", locus_tag, 
                "--outdir", output_dir, 
                "--prefix", prefix, 
                "--kingdom", cfg.config["prokka"]["kingdom"], 
                "--genus", genus,
                "--species", species,
                "--strain", strain,
                "--gcode", str(cfg.config["prokka"]["gcode"]), 
                input_file]
    return call(arguments)


def dfast_call(locus_tag, contigs_file, output_dir, sample_basename, organism):
    """
    Dfast call.
    
    Arguments:
        input_file {string} -- Contigs filename (and route).
        min_length {int} -- Minimum sequence length.
        out_path {strin} -- Output folder.
    
    Returns:
        {int} -- Execution state (0 if everything is all right)
    """
    arguments = ["dfast", 
                "--genome", contigs_file,
                "--sort_sequence", cfg.config["dfast"]["sort"],
                "--minimum_length", str(cfg.config["min_contig_len"]),
                "--use_original_name", str(cfg.config["dfast"]["use_original_name"]),
                "--step", str(cfg.config["dfast"]["step"]),
                "--organism", organism, 
                "--strain", sample_basename,
                "--locus_tag_prefix", locus_tag,
                "--out", output_dir]
    state = call(arguments)

    # Replace default output filenames including string basename
    for root, _dirs, files in os.walk(output_dir):
        for filename in files:
            if filename.__contains__("genome"):
                new_filename = filename.replace("genome", sample_basename)
                shutil.move(os.path.join(root, filename), os.path.join(root, new_filename))
    return state


def refactor_gff_from_dfast(gff_input, gff_output):
    """
    Dfast sets a generic ID in its gff file, so we replace it with one related to each sample.
    
    Arguments:
        gff_input {string} -- Input gff file (from dfast)
        gff_output {string} -- Output gff file.
    """
    with open(gff_input) as in_file, open(gff_output, 'w+') as out_file:
        for line in in_file:
            if "locus_tag=" in line:
                locus_tag = line.split("locus_tag=")[1].split(";")[0]
                line_id = line.split("\tID=")[1].split(";")[0]

                line = line.replace(line_id, locus_tag, 1)
            out_file.write(line)
        


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


def blast_postprocessing(blast_file, database_file, output_folder):
    """
    Blast post processing.
    
    Arguments:
        blast_file {string} -- Blast output file.
        database_file {string} -- Proteins database file.
        output_folder {string} -- Output folder.
    """
    
    # Read blast file as dataframe
    blast_output = pd.read_csv(blast_file, sep="\t")

    # Remove low identity rows 
    blast_output = blast_output[blast_output["% identity"] > 50]

    # Add query length and protein cover % columns
    blast_output.insert(2, "query length", "")
    blast_output.insert(3, "protein cover %", "")
    
    # Parse fasta database into dict
    proteins_database = {}
    for record in SeqIO.parse(database_file, "fasta"):
        proteins_database[record.id] = (record.description, record.seq)
    
    # Fill query length and protein cover % columns
    for index, row in blast_output.iterrows():
        query_id = row["query id"]
        query_length = len(proteins_database[query_id][1])
        blast_output.at[index, "query length"] = query_length

        # (query end - query start + 1) / query length * 100
        query_end = row["query end"]
        query_start = row["query start"]
        protein_cover = (query_end - query_start + 1) / query_length * 100
        blast_output.at[index, "protein cover %"] = protein_cover

    blast_output.to_csv(output_folder+"/BLASToutput_VF_custom_edited.txt", sep="\t",index=False)


def generate_report(samples):
    
    csv_report = pd.DataFrame(["Sample", "Reads", "AvgReadLen", "Contigs", "GenomeLen", "AvgContigLen", "N50", "GC", "DepthCov (X)"])

    for sample in samples:

        # Name of the sample
        sample = ""

        # Total number of reads after quality filtering
        n_reads = 0

        # Average read length (bp) after quality filtering
        avg_read_len = 0

        # Number of contigs of the genome (>500bp)
        n_contigs = 0

        # Length (bp) of the genome
        genome_len = 0

        # Average contig length (bp) (>500bp)
        avg_contig_len = 0

        # Length of the smallest contig in the set that contains the fewest (largest) contigs whose combined length represents at least 50% of the assembly
        n50 = 0

        # GC content (%) of the draft genome.
        gc = 0

        # Number of times each nucleotide position in the draft genome has a read that align to that position.
        depth_cov = 0


    csv_report.append({ "Sample": sample, 
                        "Reads": n_reads, 
                        "AvgReadLen": avg_read_len, 
                        "Contigs": n_contigs, 
                        "GenomeLen": genome_len, 
                        "AvgContigLen": avg_contig_len, 
                        "N50": n50, 
                        "GC": gc, 
                        "DepthCov (X)": depth_cov})


if __name__ == "__main__":

    # Get config file parameters
    annotator = cfg.config["annotator"]

    # Create output directories
    now = datetime.datetime.now()

    # Get reference files from workflow_config.py
    adapters_file =  cfg.config["adapters_reference_file"]
    reference_annotation_file = cfg.config["annotation_reference"]["file"]
    proteins_file = cfg.config["proteins_reference_file"]
    
    output_folder = sys.argv[1]

    trimmomatic_dir = output_folder+"/Trimmomatic_filtering"
    prinseq_dir = output_folder+"/Prinseq_filtering"
    flash_dir = output_folder+"/Flash_read_extension"
    spades_dir = output_folder+"/SPAdes_assembly"
    contigs_dir = output_folder+"/Contigs_renamed_shorten"
    mauve_dir = output_folder+"/Mauve_reordered_contigs"
    mauve_contigs_dir = mauve_dir+"/contigs"
    snps_dir = output_folder+"/SNP_SNIPPY"
    mlst_dir = output_folder+"/MLST"
    abricate_vir_dir = output_folder+"/ABRicate_virulence_genes"
    abricate_abr_dir = output_folder+"/ABRicate_antibiotic_resistanceGenes"
    prokka_dir = output_folder+"/Prokka_annotation"
    dfast_dir = output_folder+"/Dfast_annotation"
    roary_dir = output_folder+"/Roary_pangenome"
    roary_plots_dir = roary_dir+"/Roary_plots"
    dfast_refactor_dir = roary_dir+"/input_gff_files_edited"   # This refactor has to do with Roary so it's in Roary's folder
    blast_proteins_dir = output_folder+"/BLAST_custom_virulence_genes"
    dna_database_blast = blast_proteins_dir+"/DNA_database"

    # Create directories

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    
    if cfg.config["run_trimmomatic"]:
        os.mkdir(trimmomatic_dir)

    os.mkdir(prinseq_dir)
    os.mkdir(flash_dir)
    os.mkdir(spades_dir)
    os.mkdir(contigs_dir)
    os.mkdir(mauve_dir)
    os.mkdir(mauve_contigs_dir)
    os.mkdir(snps_dir)
    os.mkdir(mlst_dir)
    os.mkdir(abricate_vir_dir)
    os.mkdir(abricate_abr_dir)
    os.mkdir(roary_dir)

    if annotator == "dfast":
        os.mkdir(dfast_dir)
        os.mkdir(dfast_refactor_dir)
    else:
        os.mkdir(prokka_dir)

    roary_input_files = []

    # Annotate reference fasta file 
    if reference_annotation_file:
        reference_annotation_filename = reference_annotation_file.split("/")[-1]
        reference_annotation_basename = reference_annotation_filename.split(".")[-2]
        if annotator.lower() == "dfast":
            # Dfast call
            annotation_dir = dfast_dir
            print(Banner(f"\nAnnotating reference sequence: Dfast\n"), flush=True)
            dfast_call( locus_tag=reference_annotation_basename+"_L",
                        contigs_file=reference_annotation_file,
                        output_dir=dfast_dir+"/"+reference_annotation_basename,
                        sample_basename=reference_annotation_basename,
                        organism=cfg.config["annotation_reference"]["genus"]+" "+cfg.config["annotation_reference"]["species"])
            refactor_gff_from_dfast(dfast_dir+"/"+reference_annotation_basename+"/"+reference_annotation_basename+".gff",
                                    dfast_refactor_dir+"/"+reference_annotation_basename+".gff")
            # Set roary input files (renaming to get reference file first)
            os.rename(dfast_refactor_dir+"/"+reference_annotation_basename+".gff", dfast_refactor_dir+"/+"+reference_annotation_basename+".gff")
            roary_input_files.append(dfast_refactor_dir+"/+"+reference_annotation_basename+".gff")
        else:
            # Prokka call
            annotation_dir = prokka_dir
            print(Banner(f"\nAnnotating reference sequence: Prokka\n"), flush=True)
            prokka_call(locus_tag=reference_annotation_basename+"_L",
                        output_dir=prokka_dir+"/"+reference_annotation_basename,
                        prefix=reference_annotation_basename,
                        input_file=reference_annotation_file,
                        genus=cfg.config["annotation_reference"]["genus"],
                        species=cfg.config["annotation_reference"]["species"],
                        strain=cfg.config["annotation_reference"]["strain"]
                        )
            # Set roary input files (renaming to get reference file first)
            os.rename(annotation_dir+"/"+reference_annotation_basename+"/"+reference_annotation_basename+".gff",
                      annotation_dir+"/"+reference_annotation_basename+"/+"+reference_annotation_basename+".gff")
            roary_input_files.append(annotation_dir+"/"+reference_annotation_basename+"/+"+reference_annotation_basename+".gff")

    
    # Workflow Starts (for standard samples)
    sample_counter = 0
    samples_basenames = [] # Keeping track of them
    n_samples = len(read_input_files("input_files.csv"))

    for sample_basename, data in read_input_files("input_files.csv"):
        sample_fw = data["FW"]
        sample_rv = data["RV"]
        genus = data["Genus"]
        species = data["Species"]
        organism = genus+" "+species

        step_counter = 1 # Just to let the user know the number of each step
        sample_counter += 1
        samples_basenames.append(sample_basename)

        # Run trimmomatic or not
        if cfg.config["run_trimmomatic"]:
            
            print(Banner(f"\nStep {step_counter} for sequence {sample_counter}/{n_samples}({sample_basename}): Trimmomatic\n"), flush=True)
            trimmomatic_call(input_file1=sample_fw,
                            input_file2=sample_rv,
                            phred="-phred33",
                            trimfile="ILLUMINACLIP:"+adapters_file+":1:30:11",
                            paired_out_file1=trimmomatic_dir+"/"+sample_basename+"_R1_paired.fastq",
                            unpaired_out_file1=trimmomatic_dir+"/"+sample_basename+"_R1_unpaired.fastq",
                            paired_out_file2=trimmomatic_dir+"/"+sample_basename+"_R2_paired.fastq",
                            unpaired_out_file2=trimmomatic_dir+"/"+sample_basename+"_R2_unpaired.fastq")
            step_counter += 1
            prinseq_input1 = trimmomatic_dir+"/"+sample_basename+"_R1_paired.fastq"
            prinseq_input2 = trimmomatic_dir+"/"+sample_basename+"_R2_paired.fastq"
        else:
            prinseq_input1 = sample_fw
            prinseq_input2 = sample_rv


        # Create prinseq output directories
        os.mkdir(prinseq_dir+"/"+sample_basename)

        # Prinseq call
        print(Banner(f"\nStep {step_counter} for sequence {sample_counter}/{n_samples} ({sample_basename}): Prinseq\n"), flush=True)
        prinseq_call(input_file1=prinseq_input1,
                    input_file2=prinseq_input2,
                    output_folder=prinseq_dir+"/"+sample_basename,
                    sample=sample_basename,
                    log_name=prinseq_dir+"/"+sample_basename+"/"+sample_basename+".log")
        step_counter += 1
        
        # Prinseq output files refactor
        prinseq_files = refactor_prinseq_output(prinseq_dir+"/"+sample_basename, sample_basename)
        

        # Flash call
        flash_call(input_file_1=prinseq_files["R1"],
                   input_file_2=prinseq_files["R2"],
                   output_filename=sample_basename,
                   output_dir=flash_dir+"/"+sample_basename)


        # Create SPAdes output directories
        os.mkdir(spades_dir+"/"+sample_basename)

        # SPAdes call
        print(Banner(f"\nStep {step_counter} for sequence {sample_counter}/{n_samples} ({sample_basename}): SPAdes\n"), flush=True)
        
        spades_call(merged_sample=flash_dir+"/"+sample_basename+"/"+sample_basename+".extendedFrags.fastq",
                    forward_sample=flash_dir+"/"+sample_basename+"/"+sample_basename+".notCombined_1.fastq",
                    reverse_sample=flash_dir+"/"+sample_basename+"/"+sample_basename+".notCombined_2.fastq",
                    sample=sample_basename,
                    out_dir=spades_dir)
        step_counter += 1


        # Get minimum contig length
        # min_contig_threshold = get_reads_length(prinseq_dir+"/"+sample_basename+"/"+sample_basename+"_R1.fastq") * 2
        min_contig_threshold = cfg.config["min_contig_len"]

        # Trim short contigs and shorten sequences id
        contigs_trim_and_rename(contigs_file=spades_dir+"/"+sample_basename+"/"+"contigs.fasta",
                                output_filename=sample_basename+".fasta",
                                output_dir=contigs_dir,
                                min_len=min_contig_threshold)

        
        # Reordering contigs by a reference genome with MauveCM
        mauve_contigs = mauve_call(output_folder=mauve_dir+"/"+sample_basename,
                                    reference_sequence=reference_annotation_file,
                                    input_contigs=contigs_dir+"/"+sample_basename+".fasta",
                                    sample_basename=sample_basename)


        # Create Quast output directories
        quast_dir = sample_basename+"_assembly_statistics"
        os.mkdir(spades_dir+"/"+sample_basename+"/"+quast_dir)

        # Quast call
        print(Banner(f"\nStep {step_counter} for sequence {sample_counter}/{n_samples} ({sample_basename}): Quast\n"), flush=True)
        quast_call( input_file=mauve_contigs,
                    output_dir=spades_dir+"/"+sample_basename+"/"+quast_dir,
                    min_contig_len=min_contig_threshold)
        step_counter += 1


        # Annotation (Prokka or dfast)
        if annotator.lower() == "dfast":
            # Dfast call
            annotation_dir = dfast_dir
            print(Banner(f"\nStep {step_counter} for sequence {sample_counter}/{n_samples} ({sample_basename}): Dfast\n"), flush=True)
            dfast_call(locus_tag=sample_basename+"_L",
                       contigs_file=mauve_contigs,
                       output_dir=dfast_dir+"/"+sample_basename,
                       sample_basename=sample_basename,
                       organism=organism)
            step_counter += 1

            refactor_gff_from_dfast(dfast_dir+"/"+sample_basename+"/"+sample_basename+".gff", dfast_refactor_dir+"/"+sample_basename+".gff")
            # Set roary input files
            roary_input_files.append(dfast_refactor_dir+"/"+sample_basename+".gff")

        else:
            # Prokka call
            annotation_dir = prokka_dir
            print(Banner(f"\nStep {step_counter} for sequence {sample_counter}/{n_samples} ({sample_basename}): Prokka\n"), flush=True)
            prokka_call(locus_tag=sample_basename+"_L",
                        output_dir=prokka_dir+"/"+sample_basename,
                        prefix=sample_basename,
                        input_file=mauve_contigs,
                        genus=genus,
                        species=species,
                        strain=sample_basename)
            step_counter += 1
            # Set roary input files
            roary_input_files.append(annotation_dir+"/"+sample_basename+"/"+sample_basename+".gff")
        


        # SNPs identification (SNIPPY)
        print(Banner(f"\nStep {step_counter} for sequence {sample_counter}/{n_samples} ({sample_basename}): SNIPPY\n"), flush=True)
        step_counter += 1
        reference_annotation_filename = reference_annotation_file.split("/")[-1]
        reference_annotation_basename = reference_annotation_filename.split(".")[-2]
        snippy_call(reference_genome=annotation_dir+"/"+reference_annotation_basename+"/"+reference_annotation_basename+".gbk",
                    contigs=mauve_contigs,
                    output_dir=snps_dir+"/"+sample_basename,
                    prefix=sample_basename)

    # Quast report unification
    quast_report_unification(spades_dir, samples_basenames, spades_dir)

    # MLST call
    print(Banner(f"\nStep {step_counter}: MLST\n"), flush=True)
    step_counter += 1
    mlst_call(input_dir=mauve_contigs_dir,
            output_dir=mlst_dir,
            output_filename="MLST.txt")


    # ABRicate call (virulence genes)
    print(Banner(f"\nStep {step_counter}: ABRicate (virulence genes)\n"), flush=True)
    step_counter += 1
    abricate_call(input_dir=mauve_contigs_dir,
                output_dir=abricate_vir_dir,
                output_filename="VirulenceGenes.tab",
                database = cfg.config["abricate"]["virus_database"])


    # Blast call
    if cfg.config["run_blast"]:
        print(Banner(f"\nStep {step_counter}: BLAST (virulence genes against custom database)\n"), flush=True)
        step_counter += 1
        os.mkdir(blast_proteins_dir)
        os.mkdir(dna_database_blast)
        contig_files = [spades_dir+"/"+strainfolder+"/contigs.fasta" for strainfolder in next(os.walk(spades_dir))[1]]
        proteins_database_name = "VF_custom.txt"
        blast_output_name = "BLASToutput_VF_custom.txt"
        blast_call( proteins_file_ori=proteins_file, 
                    proteins_file_dest=blast_proteins_dir+"/"+proteins_database_name, 
                    contigs_files_paths=contig_files, 
                    blast_database_output=dna_database_blast+"/DNA_database.fna", 
                    blast_output_folder=blast_proteins_dir,
                    blast_output_name=blast_output_name)
        blast_postprocessing(blast_file=blast_proteins_dir+"/"+blast_output_name,
                            database_file=blast_proteins_dir+"/"+proteins_database_name,
                            output_folder=blast_proteins_dir)


    # ABRicate call (antibiotic resistance genes)
    print(Banner(f"\nStep {step_counter}: ABRicate (antibiotic resistance genes)\n"), flush=True)
    step_counter += 1
    abricate_call(input_dir=mauve_contigs_dir,
                output_dir=abricate_abr_dir,
                output_filename="AntibioticResistanceGenes.tab",
                database = cfg.config["abricate"]["bacteria_database"])
    

    # Roary call
    print(Banner(f"\nStep {step_counter}: Roary\n"), flush=True)
    step_counter += 1
    roary_call(input_files=roary_input_files, output_dir=roary_dir+"/Roary")


    # Roary plots call
    os.mkdir(roary_plots_dir)
    print(Banner(f"\nStep {step_counter}: Roary Plots\n"), flush=True)
    step_counter += 1
    roary_plots_call(input_newick=roary_dir+"/Roary/accessory_binary_genes.fa.newick",
                    input_gene_presence_absence=roary_dir+"/Roary/gene_presence_absence.csv",
                    output_dir=roary_plots_dir)

    print(Banner("\nDONE\n"), flush=True)
    # Final report

