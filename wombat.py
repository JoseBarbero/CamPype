import pandas as pd
import numpy as np
import datetime
import os
import csv
import shutil
import re
import logging
import sys
import wombat_config as cfg
import requests
from terminal_banner import Banner
from subprocess import call
from Bio import SeqIO
from io import StringIO

def welcome(wombat_img):
    img_file = open(wombat_img)
    welcome_banner =    "   _       __     __                             __       \n"\
                        "  | |     / /__  / /________  ____ ___  ___     / /_____  \n"\
                        "  | | /| / / _ \/ / ___/ __ \/ __ `__ \/ _ \   / __/ __ \ \n"\
                        "  | |/ |/ /  __/ / /__/ /_/ / / / / / /  __/  / /_/ /_/ / \n"\
                        "  |__/|__/\___/_/\___/\____/_/ /_/ /_/\___/   \__/\____/  \n"\
                        "      _       __                __          __    __\n"\
                        "     | |     / /___  ____ ___  / /_  ____ _/ /_  / /\n"\
                        "     | | /| / / __ \/ __ `__ \/ __ \/ __ `/ __/ / / \n"\
                        "     | |/ |/ / /_/ / / / / / / /_/ / /_/ / /_  /_/  \n"\
                        "     |__/|__/\____/_/ /_/ /_/_.___/\__,_/\__/ (_)   \n"    
    print(Banner(str("".join(img_file.readlines())+"\n"+welcome_banner)))

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
    arguments = ["prinseq-lite.pl", "-fastq", input_file1, "-fastq2", input_file2, "-min_len", str(cfg.config["prinseq"]["min_len"]), \
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
    arguments = ["spades.py", "--merged", merged_sample, "-1", forward_sample, "-2", reverse_sample, cfg.config["spades"]["mode"], "--cov-cutoff", cfg.config["spades"]["cov_cutoff"]]
    if cfg.config["spades"]["k"]:
        arguments.extend(["-k", str(cfg.config["spades"]["k"])])
    arguments.extend(["-o", out_dir+"/"+sample])
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
    arguments = ["MauveCM", "-output", output_folder+"/"+sample_basename, "-ref", reference_sequence, "-draft", input_contigs]    
    call(arguments)
    # Here we take the fasta file from the last iteration folder.
    shutil.copyfile(output_folder+"/"+sample_basename+"/"+max(next(os.walk(output_folder+"/"+sample_basename))[1])+"/"+sample_basename+".fasta", output_folder+"/"+sample_basename+".fasta")
    shutil.rmtree(output_folder+"/"+sample_basename)
    return output_folder+"/"+sample_basename+".fasta"


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


def mlst_call(input_dir, reference_file, output_dir, output_filename):
    """
    MLST call for every fasta file in input_dir.
    
    Arguments:
        input_dir {string} -- Input directory containing contig files.
        reference_file {string} -- Reference file directory.
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
    input_filenames.append(reference_file)
    arguments = ["mlst", *input_filenames]
    return call(arguments, stdout=output_file)


def mlst_postprocessing(mlst_file, output_file):
    col_names = ["Sample", "Genus", "ST"]
    output_data = pd.DataFrame()
    mlst_df = pd.read_csv(mlst_file, delimiter="\t",header=None)
    for _, row in mlst_df.iterrows():
        new_row = []
        new_row.append(os.path.basename(row[0]).split(".")[0])  # Basename
        new_row.append(row[1])                                  # Genus
        new_row.append(row[2])                                  # ST
        for column in mlst_df.columns[3:]:
            if not row[column].split("(")[0] in col_names:
                col_names.append(row[column].split("(")[0])
            new_row.append(int("".join(filter(str.isdigit, row[column]))))
        
        # Get clonal complex because MLST doesn't return it (https://github.com/tseemann/mlst/issues/60)
        if not "clonal_complex" in col_names:
            col_names.append("clonal_complex")
        
        url = "http://rest.pubmlst.org/db/pubmlst_campylobacter_seqdef/schemes/1/profiles_csv"
        urlData = requests.get(url).content
        database = pd.read_csv(StringIO(urlData.decode('utf-8')), sep="\t")

        new_row.append(database.loc[database["ST"] == row[2]]["clonal_complex"].values[0])


        output_data = output_data.append(pd.DataFrame([new_row], columns=col_names), ignore_index=True)

    output_data.to_csv(output_file, index=False, sep="\t")


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
    tmp_file = output_dir+"/tmp_"+output_filename
    output_file = output_dir+"/"+output_filename
    
    input_filenames = []

    for _root, _dirs, files in os.walk(input_dir):
        for filename in files:
            if filename.endswith(".fasta"):
                input_filenames.append(input_dir+"/"+filename)

    # Reference file
    input_filenames.append(cfg.config["reference_genome"]["file"])

    arguments = ["abricate", *input_filenames, "--db", database]
    
    with open(tmp_file, "w") as initial_file:
        state = call(arguments, stdout=initial_file)

    with open(tmp_file, "r") as initial_file:
        with open(output_file, "w") as final_file:
            for line in initial_file:
                if line.startswith("#FILE"):
                    line = line.replace("#FILE", "SAMPLE")
                else:
                    path = line.split("\t")[0]
                    sample = os.path.basename(path).split(".")[0]
                    line = line.replace(path, sample)
                final_file.write(line)
    os.remove(tmp_file)

    return state


def blast_call(proteins_file_ori, proteins_file_dest, contigs_files_paths, blast_database_output, blast_output_folder, blast_output_name):
    """
    Blast call.
    
    Arguments:
        proteins_file_ori {string} -- Reference proteins file path (origin).
        proteins_file_dest {string} -- Reference proteins file path (destination).
        contigs_files_paths {list} -- List of contig files.
        blast_database_output {string} -- Destination folder to blast database.
        blast_output_folder {string} -- Destination folder to blast output.
        blast_output_name {string} -- Output file name.
    
    Returns:
        {int} -- Execution state (0 if everything is all right)
    """
    # Move VF_custom.txt to ABRicateVirulenceGenes/BLAST_custom_proteins
    shutil.copy(proteins_file_ori, proteins_file_dest)

    # Concat every contig from Mauve on DNA_database.fna (replacing sequences names)
    with open(blast_database_output, "w") as output_file:
        for contig_file_path in contigs_files_paths:
            for record in SeqIO.parse(contig_file_path, "fasta"):
                strain = os.path.basename(contig_file_path).split(".")[0]
                record.id = record.name = record.description = strain+"_C_"+"_".join(record.id.split("_")[1:])
                SeqIO.write(record, output_file, "fasta")
                
    # Create blast database
    blast_db_path = os.path.dirname(os.path.abspath(blast_database_output))+"/DNA_database"
    print(blast_db_path)
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

def prokka_call(locus_tag, output_dir, prefix, input_file, genus, species, strain, proteins="", metagenome=False, rawproduct=False):
    """
    Prokka call.
    
    Arguments:
        locus_tag {string} -- Locus tag prefix.
        output_dir {string} -- Output directory.
        prefix {string} -- Filename output prefix.
        input_file {string} -- Input filename (and route).
        genus {string} -- TODO
        species {string} -- TODO
        strain {string} -- TODO
        proteins {string} -- TODO
        metagenome {string} -- TODO
        rawproduct {string} -- TODO
    
    Returns:
        {int} -- Execution state (0 if everything is all right)
    """
    arguments = ["prokka", 
                input_file,
                "--locustag", locus_tag, 
                "--outdir", output_dir, 
                "--prefix", prefix, 
                "--kingdom", cfg.config["prokka"]["kingdom"], 
                "--genus", genus,
                "--species", species,
                "--strain", strain,
                "--gcode", str(cfg.config["prokka"]["gcode"])]
    if proteins:
        arguments.extend(["--proteins", proteins])
    if metagenome:
        arguments.append("--metagenome")
    if rawproduct:
        arguments.append("--rawproduct")
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
            elif filename.__contains__("cds"):
                new_filename = filename.replace("cds", sample_basename+"_cds")
                shutil.move(os.path.join(root, filename), os.path.join(root, new_filename))
            elif filename.__contains__("protein"):
                new_filename = filename.replace("protein", sample_basename)
                shutil.move(os.path.join(root, filename), os.path.join(root, new_filename))
            elif filename.__contains__("rna"):
                new_filename = filename.replace("rna", sample_basename+"_rna")
                shutil.move(os.path.join(root, filename), os.path.join(root, new_filename))
            elif filename.__contains__("statistics"):
                new_filename = filename.replace("statistics", sample_basename+"_statistics")
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
    arguments = ["roary", "-f", output_dir, "-s", "-v", *input_files]
    ex_state = call(arguments)
    # Set Roary output directory name
    for root, dirs, _files in os.walk("."):
        for dirname in dirs:
            if dirname.startswith("Roary_pangenome_"):
                files = os.listdir(root+"/"+dirname)
                for f in files:
                    shutil.move(root+"/"+dirname+"/"+f, output_dir)
                os.rmdir(root+"/"+dirname)
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


def get_reads_table(input1, input2, sample_name, output, after_qc=False):
    """
    Gets reads info provided by prinseq-lite.
    
    Arguments:
        input1 {string} -- Forward fastq file.
        input2 {string} -- Reverse fastq file.
        sample_name {string} -- Sample name.
        output {string} -- Output file.
    
    Keyword Arguments:
        after_qc {bool} -- Describes if get_reads_table should run as if FastQC was already executed or not (default: {False})
    
    Returns:
        {dict} -- Dictionary with reads information.
    """
    arguments = ["prinseq-lite.pl", 
                "-fastq", input1, 
                "-fastq2", input2, 
                "-stats_info", 
                "-stats_len"]
    
    # Write in a bypass file prinseq-lite result
    bypass = "bypass.tsv"
    with open(bypass, "w") as f:
        state = call(arguments, stdout=f)
    call(["cat", bypass])

    # Create a new dataset with prinseq stats for every sample
    with open(bypass, "r") as bypass_file:
        data = csv.reader(bypass_file, delimiter="\t")

        has_header = False
        

        if not os.path.isfile(output):
            with open(output, "w") as out_file:
                outputwriter = csv.writer(out_file, delimiter="\t")
                if after_qc:
                    outputwriter.writerow(["Sample", "R1ReadsQC", "R1LenMeanQC", "R1LenSdQC", "R2ReadsQC", "R2LenMeanQC", "R2LenSdQC"])
                else:
                    outputwriter.writerow(["Sample", "R1Reads", "R1LenMean", "R1LenSd", "R2Reads", "R2LenMean", "R2LenSd"])
                
        
        with open(output, "a") as out_file:
            outputwriter = csv.writer(out_file, delimiter="\t")
            data_dict = {}
            
            # Uglier that indexing but more restrictive, just in case
            data_dict["Sample"] = sample_name
            for row in data:
                if row[0] == "stats_info" and row[1] == "reads":
                    data_dict["R1Reads"] = int(row[2])
                elif row[0] == "stats_info2" and row[1] == "reads":
                    data_dict["R2Reads"] = int(row[2])
                elif row[0] == "stats_len":
                    if row[1] == "stddev":
                        data_dict["R1LenSd"] = round(float(row[2]),2)
                    elif row[1] == "mean":
                        data_dict["R1LenMean"] = round(float(row[2]),2)
                elif row[0] == "stats_len2":
                    if row[1] == "stddev":
                        data_dict["R2LenSd"] = round(float(row[2]),2)
                    elif row[1] == "mean":
                        data_dict["R2LenMean"] = round(float(row[2]),2)
                        
            outputwriter.writerow([data_dict["Sample"], data_dict["R1Reads"], data_dict["R1LenMean"], data_dict["R1LenSd"], data_dict["R2Reads"], data_dict["R2LenMean"], data_dict["R2LenSd"]])
            

    call(["rm", bypass])
    return data_dict

def get_flash_reads_table(extended, notcombined1, notcombined2, sample_name, output, previous_data):
    """
    Gets reads info provided by prinseq-lite (after flash).
    
    Arguments:
        extended {string} -- Extended fastq file.
        notcombined1 {string} -- Not combined fastq forward file.
        notcombined2 {string} -- Not combined fastq reverse file.
        sample_name {string} -- Sample name.
        output {string} -- Output file.
        previous_data {dict} -- Dictionary containing reads information.
    
    Returns:
        {dict} -- Dictionary with reads information.
    """
        
    # Write in a bypass file prinseq-lite result
    bypass1 = "bypass1.tsv"
    with open(bypass1, "w") as f:
        arguments = ["prinseq-lite.pl", 
                "-fastq", extended, 
                "-stats_info", 
                "-stats_len"]
        state = call(arguments, stdout=f)
    call(["cat", bypass1])

    bypass2 = "bypass2.tsv"
    with open(bypass2, "w") as f:
        arguments = ["prinseq-lite.pl", 
                "-fastq", notcombined1, 
                "-fastq2", notcombined2, 
                "-stats_info", 
                "-stats_len"]
        state2 = call(arguments, stdout=f)
    call(["cat", bypass2])

    # Create a new dataset with prinseq stats for every sample
    with open(bypass1, "r") as bypass_file1:
        with open(bypass2, "r") as bypass_file2:
            data1 = csv.reader(bypass_file1, delimiter="\t")
            data2 = csv.reader(bypass_file2, delimiter="\t")

            has_header = False
            if not os.path.isfile(output):
                with open(output, "w") as out_file:
                    outputwriter = csv.writer(out_file, delimiter="\t")
                    outputwriter.writerow(["Sample", "JoinPercentage", "JoinReads", "JoinLenMeanReads", "JoinLenSdReads", "UnjoinR1Reads", "UnjoinR1LenMeanReads", "UnjoinR1LenSdReads", "UnjoinR2Reads", "UnjoinR2LenMeanReads", "UnjoinR2LenSdReads"])
                    
            
            with open(output, "a") as out_file:
                outputwriter = csv.writer(out_file, delimiter="\t")
                data_dict = {}
                
                data_dict["Sample"] = sample_name

                # Uglier that indexing but more restrictive, just in case
                for row in data1:
                    if row[0] == "stats_info" and row[1] == "reads":
                        data_dict["JoinReads"] = int(row[2]) * 2
                    elif row[0] == "stats_len":
                        if row[1] == "stddev":
                            data_dict["JoinLenSdReads"] = round(float(row[2]),2)
                        elif row[1] == "mean":
                            data_dict["JoinLenMeanReads"] = round(float(row[2]),2)

                data_dict["JoinPercentage"] = str(round(data_dict["JoinReads"] / (previous_data["R1Reads"] * 2) * 100, 2))+"%"

                for row in data2:
                    if row[0] == "stats_info" and row[1] == "reads":
                        data_dict["UnjoinR1Reads"] = int(row[2])
                    elif row[0] == "stats_info2" and row[1] == "reads":
                        data_dict["UnjoinR2Reads"] = int(row[2])
                    elif row[0] == "stats_len":
                        if row[1] == "stddev":
                            data_dict["UnjoinR1LenSdReads"] = round(float(row[2]),2)
                        elif row[1] == "mean":
                            data_dict["UnjoinR1LenMeanReads"] = round(float(row[2]),2)
                    elif row[0] == "stats_len2":
                        if row[1] == "stddev":
                            data_dict["UnjoinR2LenSdReads"] = round(float(row[2]),2)
                        elif row[1] == "mean":
                            data_dict["UnjoinR2LenMeanReads"] = round(float(row[2]),2)
                            
                outputwriter.writerow([ data_dict["Sample"], data_dict["JoinPercentage"], data_dict["JoinReads"], data_dict["JoinLenMeanReads"], 
                                        data_dict["JoinLenSdReads"], data_dict["UnjoinR1Reads"], data_dict["UnjoinR1LenMeanReads"], data_dict["UnjoinR1LenSdReads"], 
                                        data_dict["UnjoinR2Reads"], data_dict["UnjoinR2LenMeanReads"], data_dict["UnjoinR2LenSdReads"]])
                
    call(["rm", bypass1])
    call(["rm", bypass2])
    return data_dict


def generate_report(samples, prinseq_dir, spades_dir, mauve_dir, out_dir, info_pre_QC, info_post_QC, info_post_flash):
    """
    Creates the final report.
    
    Arguments:
        samples {list} -- Samples list.
        prinseq_dir {string} -- Prinseq results directory.
        spades_dir {string} -- SPAdes results directory.
        mauve_dir {string} --  Mauves results directory.
        out_dir {string} -- Output directory.
        infro_pre_QC {dict} -- Prinseq information before running QC.
        info_post_QC {dict} -- Prinseq information after running QC.
        info_post_flash {dict} -- Prinseq information after running flash.
    """
    
    csv_report = pd.DataFrame(columns=["Sample", "Reads", "ReadLen", "ReadsQC", "ReadsQCLen", "JoinReads", "JoinReadsLen", "Contigs", "GenomeLen", "ContigLen", "N50", "GC", "DepthCov (X)"])
    assembly_report = pd.read_csv(spades_dir+"/"+"quality_assembly_report.tsv", sep="\t")

    for sample in samples:
        
        # "Reads": Total number of reads after quality filtering
        reads = info_pre_QC[sample]["R1Reads"] + info_pre_QC[sample]["R2Reads"]
        
        # "ReadLen": Average read length (bp) before quality control.
        readlen = np.mean([info_pre_QC[sample]["R1LenMean"], info_pre_QC[sample]["R2LenMean"]])

        # "ReadsQC": Total number of reads after quality control.
        readsqc = info_post_QC[sample]["R1Reads"] + info_post_QC[sample]["R2Reads"]

        # "ReadsQCLen": Average read length (bp) after quality control.
        readsqclen = np.mean([info_post_QC[sample]["R1LenMean"], info_post_QC[sample]["R2LenMean"]])

        # "JoinReads": Total combined reads.
        joinreads = info_post_flash[sample]["JoinReads"]

        # "JoinReadsLen: Mean length of combined reads.
        joinreadslen = info_post_flash[sample]["JoinLenMeanReads"]

        # "Contigs": Number of contigs of the genome (> 500bp).
        n_contigs = int(assembly_report.loc[assembly_report['Assembly'].isin(["# contigs"])][sample])
        
        # "GenomeLen": Length (bp) of the genome.
        genome_len = int(assembly_report.loc[assembly_report['Assembly'].isin(["Total length"])][sample])
        
        # "ContigLen": Average contig length (bp) (> 500bp).
        contig_len_summatory = 0
        contig_counter = 0
        for record in SeqIO.parse(mauve_dir+"/"+sample+".fasta", "fasta"):
            contig_len_summatory += len(record.seq)
            contig_counter += 1
        avg_contig_len = contig_len_summatory/contig_counter

        # "N50": “Length of the smallest contig in the set that contains the fewest (largest) contigs whose combined length represents at least 50% of the assembly” (Miller et al., 2010).
        n50 = float(assembly_report.loc[assembly_report['Assembly'].isin(["N50"])][sample])
        
        # "GC": GC content (%) of the draft genome
        gc = float(assembly_report.loc[assembly_report['Assembly'].isin(["GC (%)"])][sample])
        
        # "DepthCov (X)": Number of times each nucleotide position in the draft genome has a read that align to that position.
        depthcov = round(info_post_flash[sample]["JoinLenMeanReads"] * info_post_flash[sample]["JoinReads"] / genome_len, 0)

        csv_report = csv_report.append({ "Sample": sample, 
                            "Reads": round(reads, 0), 
                            "ReadLen": readlen,
                            "ReadsQC": readsqc,
                            "ReadsQCLen": readsqclen, 
                            "JoinReads": joinreads, 
                            "JoinReadsLen": joinreadslen, 
                            "Contigs": n_contigs, 
                            "GenomeLen": genome_len, 
                            "ContigLen": round(avg_contig_len, 2), 
                            "N50": round(n50, 0),
                            "GC": round(gc, 2),
                            "DepthCov (X)": round(depthcov, 2)}, ignore_index=True)

    csv_report.to_csv(out_dir+"/wombat_report.csv", sep="\t", index=False)




if __name__ == "__main__":

    # Welcome
    welcome("resources/wombat_ascii.txt")

    # Get config file parameters
    annotator = cfg.config["annotator"]

    # Create output directories
    now = datetime.datetime.now()

    # Get reference files from wombat_config.py
    adapters_file =  cfg.config["adapters_reference_file"]
    reference_genome_file = cfg.config["reference_genome"]["file"]
    proteins_file = cfg.config["proteins_reference_file"]
    
    output_folder = sys.argv[1]

    trimmomatic_dir = output_folder+"/Trimmomatic_filtering"
    prinseq_dir = output_folder+"/Prinseq_filtering"
    flash_dir = output_folder+"/Flash_read_extension"
    spades_dir = output_folder+"/SPAdes_assembly"
    contigs_dir = output_folder+"/Contigs_renamed_shorten"
    mauve_dir = output_folder+"/Mauve_reordered_contigs"
    snps_dir = output_folder+"/SNP_SNIPPY"
    mlst_dir = output_folder+"/MLST"
    vir_dir = output_folder+"/Virulence_genes"
    abricate_abr_dir = output_folder+"/ABRicate_antibiotic_resistanceGenes"
    prokka_dir = output_folder+"/Prokka_annotation"
    dfast_dir = output_folder+"/Dfast_annotation"
    roary_dir = output_folder+"/Roary_pangenome"
    roary_plots_dir = roary_dir+"/Roary_plots"
    dfast_refactor_dir = roary_dir+"/input_gff_files_edited"   # This refactor has to do with Roary so it's in Roary's folder
    blast_proteins_dir = vir_dir+"/BLAST_custom_virulence_genes"
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
    os.mkdir(snps_dir)
    os.mkdir(mlst_dir)
    os.mkdir(vir_dir)
    os.mkdir(abricate_abr_dir)
    os.mkdir(roary_dir)

    if annotator == "dfast":
        os.mkdir(dfast_dir)
        os.mkdir(dfast_refactor_dir)
    else:
        os.mkdir(prokka_dir)

    roary_input_files = []
    summary_pre_qc = {}
    summary_post_qc = {}
    summary_post_flash = {}

    # Annotate reference fasta file 
    if reference_genome_file:
        reference_genome_filename = reference_genome_file.split("/")[-1]
        reference_genome_basename = reference_genome_filename.split(".")[-2]
        if annotator.lower() == "dfast":
            # Dfast call
            annotation_dir = dfast_dir
            print(Banner(f"\nAnnotating reference sequence: Dfast\n"), flush=True)
            dfast_call( locus_tag=reference_genome_file+"_L",
                        contigs_file=reference_genome_file,
                        output_dir=dfast_dir+"/"+reference_genome_basename,
                        sample_basename=reference_genome_basename,
                        organism=cfg.config["reference_genome"]["genus"]+" "+cfg.config["reference_genome"]["species"])
            refactor_gff_from_dfast(dfast_dir+"/"+reference_genome_basename+"/"+reference_genome_basename+".gff",
                                    dfast_refactor_dir+"/"+reference_genome_basename+".gff")
            # Set roary input files (renaming to get reference file first)
            os.rename(dfast_refactor_dir+"/"+reference_genome_basename+".gff", dfast_refactor_dir+"/+"+reference_genome_basename+".gff")
            roary_input_files.append(dfast_refactor_dir+"/+"+reference_genome_basename+".gff")
        else:
            # Prokka call
            annotation_dir = prokka_dir
            print(Banner(f"\nAnnotating reference sequence: Prokka\n"), flush=True)
            prokka_call(locus_tag=reference_genome_basename+"_L",
                        output_dir=prokka_dir+"/"+reference_genome_basename,
                        prefix=reference_genome_basename,
                        input_file=reference_genome_file,
                        genus=cfg.config["reference_genome"]["genus"],
                        species=cfg.config["reference_genome"]["species"],
                        strain=cfg.config["reference_genome"]["strain"],
                        proteins=cfg.config["reference_genome"]["proteins"],
                        metagenome=False, # False in reference file
                        rawproduct=cfg.config["prokka"]["rawproduct"]
                        )
            # Set roary input files (renaming to get reference file first)
            os.rename(annotation_dir+"/"+reference_genome_basename+"/"+reference_genome_basename+".gff",
                      annotation_dir+"/"+reference_genome_basename+"/+"+reference_genome_basename+".gff")
            roary_input_files.append(annotation_dir+"/"+reference_genome_basename+"/+"+reference_genome_basename+".gff")

    
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
            
            print(Banner(f"\nStep {step_counter} for sequence {sample_counter}/{n_samples} ({sample_basename}): Trimmomatic\n"), flush=True)
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
        print(Banner(f"\nStep {step_counter} for sequence {sample_counter}/{n_samples} ({sample_basename}): Flash\n"), flush=True)
        flash_call(input_file_1=prinseq_files["R1"],
                   input_file_2=prinseq_files["R2"],
                   output_filename=sample_basename,
                   output_dir=flash_dir+"/"+sample_basename)
        step_counter += 1

        # Quality reports
        print(Banner(f"\nStep {step_counter} for sequence {sample_counter}/{n_samples} ({sample_basename}): Read statistics\n"), flush=True)
        step_counter += 1
        report_pre_qc = get_reads_table(sample_fw, sample_rv, sample_basename, prinseq_dir+"/reads_statistics_beforeQC.tsv", False)

        report_post_qc = get_reads_table(prinseq_files["R1"], prinseq_files["R2"], sample_basename, prinseq_dir+"/reads_statistics_afterQC.tsv", True)

        report_post_flash = get_flash_reads_table(flash_dir+"/"+sample_basename+"/"+sample_basename+".extendedFrags.fastq", 
                                flash_dir+"/"+sample_basename+"/"+sample_basename+".notCombined_1.fastq",
                                flash_dir+"/"+sample_basename+"/"+sample_basename+".notCombined_2.fastq",
                                sample_basename, flash_dir+"/reads_statistics_FLASH.tsv", report_post_qc)

        summary_pre_qc[sample_basename] = report_pre_qc
        summary_post_qc[sample_basename] = report_post_qc
        summary_post_flash[sample_basename] = report_post_flash

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
        print(Banner(f"\nStep {step_counter} for sequence {sample_counter}/{n_samples} ({sample_basename}): reordering {sample_basename} genome against reference genome\n"), flush=True)
        mauve_contigs = mauve_call(output_folder=mauve_dir,
                                    reference_sequence=reference_genome_file,
                                    input_contigs=contigs_dir+"/"+sample_basename+".fasta",
                                    sample_basename=sample_basename)
        step_counter += 1

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
                        strain=sample_basename,
                        proteins=cfg.config["reference_genome"]["proteins"],
                        metagenome=cfg.config["prokka"]["metagenome"],
                        rawproduct=cfg.config["prokka"]["rawproduct"])
            step_counter += 1
            # Set roary input files
            roary_input_files.append(annotation_dir+"/"+sample_basename+"/"+sample_basename+".gff")
        


        # SNPs identification (SNIPPY)
        print(Banner(f"\nStep {step_counter} for sequence {sample_counter}/{n_samples} ({sample_basename}): SNIPPY\n"), flush=True)
        step_counter += 1
        reference_genome_filename = reference_genome_file.split("/")[-1]
        reference_genome_basename = reference_genome_filename.split(".")[-2]
        snippy_call(reference_genome=cfg.config["reference_genome"]["proteins"],
                    contigs=mauve_contigs,
                    output_dir=snps_dir+"/"+sample_basename,
                    prefix=sample_basename)

    # Quast report unification
    quast_report_unification(spades_dir, samples_basenames, spades_dir)

    # MLST call
    mlst_out_file = "MLST.txt"
    print(Banner(f"\nStep {step_counter}: MLST\n"), flush=True)
    step_counter += 1
    mlst_call(input_dir=mauve_dir,
            reference_file=reference_genome_file,
            output_dir=mlst_dir,
            output_filename=mlst_out_file)
    
    # MLST postprocessing
    mlst_postprocessing(mlst_dir+"/"+mlst_out_file, mlst_dir+"/MLST_edited.txt")

    # ABRicate call (virulence genes)
    print(Banner(f"\nStep {step_counter}: ABRicate (virulence genes)\n"), flush=True)
    step_counter += 1
    abricate_call(input_dir=mauve_dir,
                output_dir=vir_dir,
                output_filename="VirulenceGenes.tab",
                database = cfg.config["abricate"]["virus_database"])


    # Blast call
    if cfg.config["run_blast"]:
        print(Banner(f"\nStep {step_counter}: BLAST (virulence genes against custom database)\n"), flush=True)
        step_counter += 1
        contigs_dir = mauve_dir
        contig_files = ([os.path.join(contigs_dir, f) for f in os.listdir(contigs_dir)])
        contig_files.append(cfg.config["reference_genome"]["file"])

        blast_output_name = "BLASToutput_VF_custom.txt"
        proteins_file = cfg.config["proteins_reference_file"]
        dna_database_blast = blast_proteins_dir+"/DNA_database"
        proteins_database_name = "VF_custom.txt"    # This is an output file name
        if not os.path.exists(blast_proteins_dir):
            os.mkdir(blast_proteins_dir)
        if not os.path.exists(dna_database_blast):
            os.mkdir(dna_database_blast)

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
    abricate_call(input_dir=mauve_dir,
                output_dir=abricate_abr_dir,
                output_filename="AntibioticResistanceGenes.tab",
                database = cfg.config["abricate"]["bacteria_database"])
    

    # Roary call
    print(Banner(f"\nStep {step_counter}: Roary\n"), flush=True)
    step_counter += 1
    roary_call(input_files=roary_input_files, output_dir=roary_dir)


    # Roary plots call
    os.mkdir(roary_plots_dir)
    print(Banner(f"\nStep {step_counter}: Roary Plots\n"), flush=True)
    step_counter += 1
    roary_plots_call(input_newick=roary_dir+"/accessory_binary_genes.fa.newick",
                    input_gene_presence_absence=roary_dir+"/gene_presence_absence.csv",
                    output_dir=roary_plots_dir)

    print(Banner("\nDONE\n"), flush=True)
    
    
    # Final report
    generate_report(samples_basenames, 
                    prinseq_dir, 
                    spades_dir, 
                    mauve_dir, 
                    output_folder, 
                    summary_pre_qc, summary_post_qc, summary_post_flash)
