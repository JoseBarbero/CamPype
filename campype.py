import numpy as np
import pandas as pd
import datetime
import os
import csv
import shutil
import re
import logging
import copy
import sys
import campype_config as cfg
import requests
import json
from terminal_banner import Banner
from subprocess import call
from Bio import SeqIO
from io import StringIO

def welcome(banner_img):
    img_file = open(banner_img)
    print(Banner(str("".join(img_file.readlines()))), flush=True)

def read_input_files(indexfile):
    """
    Gets every pair of reads on input_files.csv
    
    Arguments:
        indexfile {string} -- Filename (and route) of the file containing inputs.
    
    Returns:
        files_data {list of tuples} -- List of tuples each containing forward read file path, reverse read file path and file basename (just the sample number).
    """
    files_df = pd.read_csv(indexfile, sep="\t")
    # Replace nans with empty strings
    files_df = files_df.fillna("")
    files_data = []
    for _, row in files_df.iterrows():
        files_data.append((str(row["Samples"]), {"FW": str(row["Forward"]), 
                                                 "RV": str(row["Reverse"]), 
                                                 "Genus": str(row["Genus"]), 
                                                 "Species": str(row["Species"])}))
    return files_data

def trimmomatic_call(input_file1, input_file2, phred, trimfile,
                    paired_out_file1, paired_out_file2, unpaired_out_file1, unpaired_out_file2, threads=None):
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
    if threads:
        arguments.extend(["-threads", str(threads)])
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
    arguments = ["prinseq-lite.pl", "-fastq", input_file1, "-fastq2", input_file2, "-min_len", str(cfg.config["read_qc_filtering"]["min_len"]), \
                "-min_qual_mean", str(cfg.config["read_qc_filtering"]["min_qual_mean"]), "-trim_qual_right", str(cfg.config["read_qc_filtering"]["trim_qual_right"]), "-trim_qual_window", \
                str(cfg.config["read_qc_filtering"]["trim_qual_window"]), "-trim_qual_type", "mean", "-out_format", "3", "-out_good", output_folder+"/"+sample, "-out_bad", "null", log_name]
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
    
    # If we are in compressed mode we need to compress the outputs
    if compressed_mode:
        for key, value in filenames.items():
            arguments = ["gzip", value, "-c"]
            call(arguments, stdout=open(value+".gz", 'w'))
    return filenames

def kraken_call(db_file, output_file, fw_file, rv_file=None, threads=None):
    """
    Kraken call.

    Arguments:
        db_file {string} -- Kraken database file.
        output_file {string} -- Output file (report).
        fw_file {string} -- Forward file (from Prinseq if we are working in fastq mode)
        rv_file {string} -- Reverse prinseq file (from Prinseq if we are working in fastq mode)
    """

    # If we don't have reverse reads then we are working with the assembled genome
    if rv_file:
        arguments = ["kraken2", "--db", db_file, "--paired", fw_file, rv_file, "--output", "-", "--report", output_file]
    else:
        arguments = ["kraken2", "--db", db_file, fw_file, "--output", "-", "--report", output_file]

    if threads:
        arguments.extend(["--threads", str(threads)])

    call(arguments)

    # Get genus and species from kraken output file
    kraken_report_columns = ["Percentage", "Reads_clade_covered", "Reads_clade_assigned", "TaxRank", "NCBI_taxID", "Scientific_name"]
    sample_report_df = pd.read_csv(output_file, sep="\t", header=None, names=kraken_report_columns)
    # Only keep rows with species information
    sample_report_df = sample_report_df[sample_report_df["TaxRank"] == "S"]
    # Keep only the first row (the one with the highest percentage)
    sample_report_df = sample_report_df.head(1)

    genus = sample_report_df["Scientific_name"].values[0].split(" ")[0]
    species = sample_report_df["Scientific_name"].values[0].split(" ")[1]

    return genus, species

def kraken_report_unification(kraken_reports, output_file):
    """
    Creates a report unifying reports from every sample.
    
    Arguments:
        kraken_reports {list} -- List of tuples of kraken reports (sample, report file).
        output_file {string} -- Output file.
    """
    kraken_report_columns = ["Percentage", "Reads_clade_covered", "Reads_clade_assigned", "TaxRank", "NCBI_taxID", "Scientific_name"]
    combined_report_df = pd.DataFrame(columns=["Sample", "Species", "Species_percentage"])
    for sample, sample_report in kraken_reports:
        sample_report_df = pd.read_csv(sample_report, sep="\t", header=None, names=kraken_report_columns)
        # Only keep rows with species information
        sample_report_df = sample_report_df[sample_report_df["TaxRank"] == "S"]
        # Keep only the first row (the one with the highest percentage)
        sample_report_df = sample_report_df.head(1)
        # Add to combined report
        scientific_name = sample_report_df["Scientific_name"].values[0].strip() # Kraken uses spaces to represent the hierarchy so we need to strip them
        species_percentage = sample_report_df["Percentage"].values[0]
        combined_report_df = combined_report_df.append({"Sample": sample, "Species": scientific_name, "Species_percentage": species_percentage}, ignore_index=True)
    
    # Export report to tsv
    combined_report_df.to_csv(output_file, sep="\t", index=False, header=True)
    
    
def flash_call(input_file_1, input_file_2, output_filename, output_dir, threads=None):
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

    if threads:
        arguments.extend(["--threads", str(threads)])
    return call(arguments)


def spades_call(forward_sample, reverse_sample, sample, out_dir, merged_sample=None, threads=None):
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
    if merged_sample:
        arguments = ["spades.py", "--merged", merged_sample, "-1", forward_sample, "-2", reverse_sample, cfg.config["assembly"]["mode"]]
    else:
        arguments = ["spades.py", "-1", forward_sample, "-2", reverse_sample, cfg.config["assembly"]["mode"]]

    if cfg.config["assembly"]["k"]:
        arguments.extend(["-k", str(cfg.config["assembly"]["k"])])
    arguments.extend(["-o", out_dir+"/"+sample])

    if threads:
        arguments.extend(["--threads", str(threads)])

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
    contigs_basename = os.path.basename(input_contigs).split(".")[0]   # Mauve names it's output files after the input file basename
    fasta_output = max(os.walk(output_folder+"/"+sample_basename))[0]+"/"+[f for f in max(os.walk(output_folder+"/"+sample_basename))[2] if f.startswith(contigs_basename) and (f.endswith('.fna') or f.endswith('.fasta'))][0]
    shutil.copyfile(fasta_output, output_folder+"/"+sample_basename+".fasta")
    shutil.rmtree(output_folder+"/"+sample_basename)

    # Mauve sets locus string too long to process to tools like dfast. So we need to reprocess those loci to make them shorter.
    shorten_loci(output_folder+"/"+sample_basename+".fasta")

    return output_folder+"/"+sample_basename+".fasta"

def shorten_loci(original_file):
    """
    For BioPython loci must have names lower than 16 characters.
    This function will remove everything after the second "_" character.

    Args:
        originalfile (str): Path to original file.
        correctedfile (str): Path to corrected file.
    """
    corrected_file = original_file+".tmp"
    with open(original_file) as original, open(corrected_file, 'w') as corrected:
        records = SeqIO.parse(original_file, 'fasta')
        for record in records:
            original_id = record.id
            record.id = "_".join(original_id.split("_")[:2])
            record.description = "_".join(original_id.split("_")[2:])
            SeqIO.write(record, corrected, 'fasta')
    os.remove(original_file)
    shutil.move(corrected_file, original_file)

def snippy_call(reference_genome, contigs, output_dir, prefix, threads=None):
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

    if threads:
        arguments.extend(["--cpus", str(threads)])

    return call(arguments)


def snippy_summary(snippy_files, output_file):
    """
    Generates a summary table from snippy data.

    Args:
        snippy_files {list}: List of snippy files.
        output_file {str}: Summary output file.
    """
    summary_df = pd.DataFrame(columns=["Sample", "Variant-COMPLEX", "Variant-DEL", "Variant-INS", "Variant-SNP", "VariantTotal"])
    for snippy_file in snippy_files:
        snippy_data = {}

        # Get snippy data into a dictionary
        with open(snippy_file) as f:
            for line in f:
                snippy_data[line.split("\t")[0]] = line.split("\t")[1].replace("\n", "")
        
        new_row = {}
        new_row["Sample"] = os.path.splitext(os.path.basename(snippy_file))[0]  # We get the file basename to get the sample name
        if 'Variant-COMPLEX' in snippy_data:
            new_row["Variant-COMPLEX"] = snippy_data["Variant-COMPLEX"]
        if 'Variant-DEL' in snippy_data:
            new_row["Variant-DEL"] = snippy_data["Variant-DEL"]
        if 'Variant-INS' in snippy_data:
            new_row["Variant-INS"] = snippy_data["Variant-INS"]
        if 'Variant-SNP' in snippy_data:
            new_row["Variant-SNP"] = snippy_data["Variant-SNP"]
        if 'VariantTotal' in snippy_data:
            new_row["VariantTotal"] = snippy_data["VariantTotal"]

        summary_df = summary_df.append(new_row, ignore_index=True)

        # Export to tsv
        summary_df.to_csv(output_file, sep="\t",index=False)


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
        if len(record.seq) >= min_len:
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


def quast_call(input_file, output_dir, min_contig_len, threads=None):
    """
    Quast call.
    
    Arguments:
        input_file {string} -- Input file (and route).
        output_dir {string} -- Output directory.
        min_contig_len -- Lower threshold for a contig length (in bp).
    
    Returns:
        {int} -- Execution state (0 if everything is all right)
    """
    arguments = ["quast", input_file, "-o", output_dir, "--min-contig", str(min_contig_len), "--no-icarus", "--silent"]

    if threads:
        arguments.extend(["--threads", str(threads)])

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


def mlst_call(input_dir, reference_file, output_dir, output_filename, threads=None):
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

    if threads:
        arguments.extend(["--threads", str(threads)])

    return call(arguments, stdout=output_file)


def mlst_postprocessing(mlst_file, output_file):
    col_names = ["Sample", "Genus", "ST"]
    output_data = pd.DataFrame()
    mlst_df = pd.read_csv(mlst_file, delimiter="\t", header=None)

    url = "http://rest.pubmlst.org/db/pubmlst_campylobacter_seqdef/schemes/1/profiles_csv"
    urlData = requests.get(url).content
    database = pd.read_csv(StringIO(urlData.decode('utf-8')), sep="\t")

    for _, row in mlst_df.iterrows():
        new_row = []
        new_row.append(os.path.basename(row[0]).split(".fast")[0])  # Basename
        new_row.append(row[1])                                  # Genus
        new_row.append(row[2])                                  # ST
        for column in mlst_df.columns[3:]:
            if not row[column].split("(")[0] in col_names:
                col_names.append(row[column].split("(")[0])
            new_row.append(int("".join(filter(str.isdigit, row[column]))))
        
        # Get clonal complex because it's not returned by MLST (https://github.com/tseemann/mlst/issues/60)
        if not "clonal_complex" in col_names:
            col_names.append("clonal_complex")
        
        try:
            st = int(row[2])
            if len(database.loc[database["ST"] == st]["clonal_complex"].values) > 0:
                new_row.append(database.loc[database["ST"] == st]["clonal_complex"].values[0])
            else:
                new_row.append("Other")
        except ValueError:
            new_row.append("Other")
        


        output_data = output_data.append(pd.DataFrame([new_row], columns=col_names), ignore_index=True)
        # Sample column is a string
        output_data["Sample"] = output_data["Sample"].astype(str)
        print('MLST post processing:')
        print(output_data)
    output_data.to_csv(output_file, index=False, sep="\t")


def abricate_call(input_dir, output_dir, output_filename, database, mincov=False, minid=False, gene_matrix_file=False, samples=False, threads=None):
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
    output_matrix = output_dir+"/"+gene_matrix_file
    
    input_filenames = []

    for _root, _dirs, files in os.walk(input_dir):
        for filename in files:
            if filename.endswith(".fasta"):
                input_filenames.append(input_dir+"/"+filename)

    arguments = ["abricate", *input_filenames, "--db", database]
    
    if mincov:
        arguments.extend(["--mincov", str(mincov)])
    if minid:
        arguments.extend(["--minid", str(minid)])

    if threads:
        arguments.extend(["--threads", str(threads)])
    
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
    
    if gene_matrix_file:
        abricate_presence_absence_matrix(output_file, output_matrix, samples, database)

    return state


def abricate_presence_absence_matrix(abricate_file, p_a_matrix_file, samples, database):
    """
    Creates the presence/absence matrix for the input genes file.

    Arguments:
        abricate_file {string} -- File from abricate.
        p_a_matrix_file {string} -- Gene presence/absence matrix.
        samples {list} -- List of samples to analyze.
        database {str} -- Database.
    """

    abricate_data = pd.read_csv(abricate_file, sep="\t")
    abricate_data["SAMPLE"] = abricate_data["SAMPLE"].astype(str)    # SAMPLE column in abricate_data is a string
    gene_presence_absence = pd.DataFrame(columns=["Gene", *samples, "RESISTANCE", "DATABASE"])

    genes_type_ABR = dict(zip(abricate_data["GENE"], abricate_data["RESISTANCE"]))
    genes_type_VIR = dict(zip(abricate_data["GENE"], abricate_data["PRODUCT"]))
    for gene in genes_type_ABR.keys():
        gene_content = abricate_data[abricate_data["GENE"] == gene]
        new_row = {"Gene": gene}
        new_row["RESISTANCE"] = str(genes_type_ABR[gene]).lower()
        if new_row["RESISTANCE"] == "nan":
            new_row["RESISTANCE"] = str(genes_type_VIR[gene]).lower()
        for sample in samples:
            sample_content = gene_content[gene_content["SAMPLE"] == sample]
            if len(sample_content) == 1:
                if sample_content["%COVERAGE"].max() >= cfg.config["virulence_genes"]["abricate"]["mincov"] and sample_content["%IDENTITY"].max() >= cfg.config["virulence_genes"]["abricate"]["minid"]:
                    new_row[sample] = 1
                else:
                    new_row[sample] = 0
            else:
                new_row[sample] = 0
        new_row["DATABASE"] = database
        gene_presence_absence = gene_presence_absence.append(new_row, ignore_index=True)

    # Sort rows
    gene_presence_absence = gene_presence_absence.sort_values(by=["RESISTANCE", "Gene"])
    
    # Export to tsv
    gene_presence_absence.to_csv(p_a_matrix_file, sep="\t",index=False)



def amrfinder_call(amrfinder_out_file, resume_file, samples_basenames, annotation_dir, gff_dir, genus, db_name, output_dir, ref_genome_basename=None, threads=None):
    has_header = False

    samples = samples_basenames

    with open(amrfinder_out_file, "w") as global_file:
        for sample in samples:
            
            gff_file = sample
            # Ref genome sample starts by + so that roary can place it first
            if sample == ref_genome_basename:             
                gff_file = "+"+sample
            arguments = ["amrfinder",
                        "-n", annotation_dir+"/"+sample+"/"+sample+".fna", 
                        "-p", annotation_dir+"/"+sample+"/"+sample+".faa", 
                        "-g", gff_dir+"/"+gff_file+".gff",
                        "--plus", 
                        "-i", str(cfg.config["antimicrobial_resistance_genes"]["amrfinder"]["minid"]/100),
                        "-c", str(cfg.config["antimicrobial_resistance_genes"]["amrfinder"]["mincov"]/100),
                        "-o", output_dir+"/"+sample+".txt"]
            if genus == "Campylobacter":
                arguments.extend(["--organism", genus])
            if threads:
                arguments.extend(["--threads", str(threads)])
                
            call(arguments)

            # Group all the results in a single file        
            with open(output_dir+"/"+sample+".txt") as in_file:
                lines = in_file.readlines()
                if has_header == False:
                    global_file.write("Sample\t"+lines[0])
                    has_header = True
                for line in lines[1:]:
                    global_file.write(sample+"\t"+line)
            os.remove(output_dir+"/"+sample+".txt")
    
    # Generate AMR Genes Resume
    amr_data = pd.read_csv(amrfinder_out_file, sep="\t")
    columns = amr_data["Class"].unique().tolist()
    columns.append("DATABASE")
    
    resume = pd.DataFrame(columns=columns)
    samples_data = {}
    for sample in samples:
        samples_data[sample] = {"Sample": sample}
    for _, row in amr_data[(amr_data["Element type"] == "AMR") & (amr_data["Element subtype"] == "AMR")].iterrows():
        for column in columns:
            if row["Class"] == column:
                samples_data[str(row["Sample"])][column] = row["Gene symbol"]
            
    for sample, values in samples_data.items():
        resume = resume.append(values, ignore_index=True)
    
    # Sort columns
    resume = resume[["Sample"] + columns]
    resume = resume.dropna(how ="all", axis=1)
    resume = resume.fillna("-")

    resume["DATABASE"] = db_name

    resume.to_csv(resume_file, sep="\t", index=False)
    with open(resume_file, "a") as res_file:
        res_file.write("'-' means no AMR gene found")

def amrfinder_get_point_mutations(amrfinder_file, output_file):
    """
    Processes the AMR Finder matrix to get the point mutations.
    
    We get only the results from the column "Element subtype" when the value is "POINT".
    In this new table, the first column will be Sample and will be filled with the names of the different Samples.
    The headers of the following columns will be the different contents of the column "Class" (when "Element subtype" = "POINT").
    The table will be filled with the content of the column "Gene symbol" when the Sample matches Class. BUT, a same Class can have more than one Gene symbol different.
    """

    # Read the matrix 
    amr_data = pd.read_csv(amrfinder_file, sep="\t")

    # Get the point mutations lines
    point_mutations = amr_data[amr_data["Element subtype"] == "POINT"]

    # Get the different classes
    classes = point_mutations["Class"].unique().tolist()

    # Set the header of the new table
    header = ["Sample"] + classes
    
    # Create the new table
    point_mutations_table = pd.DataFrame(columns=header)

    # Get the samples
    samples = amr_data["Sample"].unique().tolist()

    for sample in samples:
        point_mutations_table.loc[sample] = sample
        for sample_class in classes:
            gene_symbols = point_mutations[(point_mutations["Sample"] == sample) & (point_mutations["Class"] == sample_class)]["Gene symbol"].unique().tolist()
            if len(gene_symbols) == 1:
                point_mutations_table.loc[sample, sample_class] = gene_symbols[0]
            elif len(gene_symbols) > 1:
                point_mutations_table.loc[sample, sample_class] = ", ".join(gene_symbols)
            else:
                point_mutations_table.loc[sample, sample_class] = "-"

    # Export the new matrix
    point_mutations_table.to_csv(output_file, sep="\t", index=False)

    # Add footer
    with open(output_file, "a") as res_file:
        res_file.write("'-' means no point mutation found")
    

def blast_call(proteins_file_ori, proteins_file_dest, contigs_files_paths, blast_database_output, blast_output_folder, blast_output_name, threads=None):
    """
    Blast call.
    
    Arguments:
        proteins_file_ori {string} -- Reference proteins file path (origin).
        proteins_file_dest {string} -- Reference proteins file path (destination).create new partition windows 10
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
    call(["makeblastdb", "-in", blast_database_output, "-dbtype", "nucl", 
          "-out", blast_db_path, "-title", "DNA_Database"])

    # Call tblastn
    arguments = ["tblastn", "-db", blast_db_path, "-query", proteins_file_dest,
                "-soft_masking", str(cfg.config["virulence_genes"]["blast"]["soft_masking"]).lower(),
                "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq",
                "-out", blast_output_folder+"/"+blast_output_name]
    if threads:
        arguments.extend(["-num_threads", str(threads)])

    tblastn_state = call(arguments)

    # Add header to tblastn output
    with open(blast_output_folder+"/"+blast_output_name, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        headers = ["query id", "subject id", "% identity", "alignment length", "mismatches", "gap opens", 
                   "query start", "query end", "subject start", "subject end", "evalue", "score", "aligned part of subject sequence"]

        f.write("\t".join(headers) + '\n' + content)

    return tblastn_state


def blast_postprocessing(blast_file, database_file, output_folder, samples):
    """
    Blast post processing.
    
    Arguments:
        blast_file {string} -- Blast output file.
        database_file {string} -- Proteins database file.
        output_folder {string} -- Output folder.
        samples {[string]} -- List of samples names.
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

    # Split subject_id column
    blast_output[["Sample", "Contig sample"]] = blast_output["subject id"].str.split("_", 1, expand=True)
    blast_output = blast_output.drop("subject id", axis=1)
    
    # Rename columns
    blast_output.rename(columns={"query id": "Protein",
                                "query length": "Protein length",
                                "protein cover %": "% protein cover", 
                                "% identity": "% protein identity",
                                "alignment length": "Alignment length",
                                "mismatches": "Mismatches",
                                "gap opens": "Gap opens",
                                "query start": "Protein start",
                                "query end": "Protein end",
                                "subject start": "Sample start",
                                "subject end": "Sample end",
                                "evalue": "Evalue",
                                "score": "Score",
                                "aligned part of subject sequence": "Aligned part of subject sequence"}, 
                                inplace=True)

    # Add protein type column
    proteins_dict = {}
    for record in SeqIO.parse(database_file, "fasta"):
        proteins_dict[record.id] = record.description.split()[1]
    blast_output["Protein type"] = blast_output.apply (lambda row: proteins_dict[row["Protein"]].lower(), axis=1)

    # Sort columns
    new_columns = [ "Protein",
                    "Protein type",
                    "Sample",
                    "Contig sample",
                    "Protein length",
                    "% protein cover",
                    "% protein identity",
                    "Alignment length",
                    "Mismatches",
                    "Gap opens",
                    "Protein start",
                    "Protein end",
                    "Sample start",
                    "Sample end",
                    "Evalue",
                    "Score",
                    "Aligned part of subject sequence"]
    
    blast_output = blast_output[new_columns]
    blast_output['% protein cover'] = blast_output['% protein cover'].apply(lambda x: round(x, 2))
    blast_output['% protein identity'] = blast_output['% protein identity'].apply(lambda x: round(x, 2))

    # Sort rows
    blast_output = blast_output.sort_values(by=["Protein type", "Protein", "Sample"])
    # Sample column must be a string
    blast_output["Sample"] = blast_output["Sample"].astype(str)

    # Export to tsv
    blast_output.to_csv(output_folder+"/Virulence_genes_BLAST_processed.tsv", sep="\t",index=False)    

    # Create presence/absence matrix
    get_presence_absence_matrix(samples, proteins_dict, blast_output, output_folder+"/Virulence_genes_BLAST_matrix.tsv")


def get_presence_absence_matrix(samples, genes_type, blast_df, p_a_matrix_file):
    """
    Creates the presence/absence matrix for the input genes file.

    Arguments:
        samples {[strings]} -- Samples names list.
        genes_type {dict} -- Genes and its type dict.
        blast_df {pandas dataframe} -- Sample/genes information from blast.
        p_a_matrix_file {string} -- Gene presence/absence matrix.
    """

    gene_presence_absence = pd.DataFrame(columns=["Protein", *samples, "Type"])

    for gene in genes_type.keys():
        gene_content = blast_df[blast_df["Protein"] == gene]
        new_row = {"Protein": gene}
        new_row["Type"] = genes_type[gene].lower()
        for sample in samples:
            sample_content = gene_content[gene_content["Sample"] == sample]
            if len(sample_content) > 0:
                if sample_content["% protein cover"].max() >= cfg.config["virulence_genes"]["blast"]["presence_absence_matrix"]["mincov"] and sample_content["% protein identity"].max() >= cfg.config["virulence_genes"]["blast"]["presence_absence_matrix"]["minid"]:
                    new_row[sample] = 1
                else:
                    new_row[sample] = 0
            else:
                new_row[sample] = 0
        gene_presence_absence = gene_presence_absence.append(new_row, ignore_index=True)
    
    # capA special distinction (the protein is split in two fragments)
    capA_row = {"Protein": "capA", "Type": genes_type["capA_ORF1"].lower()}
    orf1_df = gene_presence_absence[gene_presence_absence["Protein"] == "capA_ORF1"]
    orf2_df = gene_presence_absence[gene_presence_absence["Protein"] == "capA_ORF2"]
    for sample in samples:
        if orf1_df.iloc[0][sample] == 1 and orf2_df.iloc[0][sample] == 1:
            capA_row[sample] = 1
        else:
            capA_row[sample] = 0

    gene_presence_absence = gene_presence_absence.append(capA_row, ignore_index=True)
    gene_presence_absence = gene_presence_absence[gene_presence_absence["Protein"] != "capA_ORF1"]
    gene_presence_absence = gene_presence_absence[gene_presence_absence["Protein"] != "capA_ORF2"]
    gene_presence_absence["DATABASE"] = "Campylobacter_custom_VFDB"

    # Sort rows
    gene_presence_absence = gene_presence_absence.sort_values(by=["Type", "Protein"])
    
    # Export to tsv
    gene_presence_absence.to_csv(p_a_matrix_file, sep="\t",index=False)
    with open(p_a_matrix_file, "a") as matrix_file:
        matrix_file.write("Coverage >= " + 
                            str(cfg.config["virulence_genes"]["blast"]["presence_absence_matrix"]["mincov"]) +
                            " % and identity >= " +
                            str(cfg.config["virulence_genes"]["blast"]["presence_absence_matrix"]["minid"]) + 
                            " % on each sample for considering a virulence gene as present.")

def prokka_call(locus_tag, output_dir, prefix, input_file, genus, species, strain, proteins="", rawproduct=False, threads=None):
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
        rawproduct {string} -- TODO
    
    Returns:
        {int} -- Execution state (0 if everything is all right)
    """
    arguments = ["prokka", 
                input_file,
                "--locustag", locus_tag, 
                "--outdir", output_dir, 
                "--prefix", prefix, 
                "--kingdom", "Bacteria", 
                "--genus", genus,
                "--species", species,
                "--strain", strain,
                "--gcode", "11"]
    if cfg.config["annotation"]["prokka"]["reference_annotation"]:
        arguments.extend(["--proteins", proteins])
    if rawproduct:
        arguments.append("--rawproduct")
    if threads:
        arguments.extend(["--cpus", str(threads)])

    return call(arguments)


def prokka_summary(prokka_files, output_file):
    """
    Generates a summary table from prokka data.

    Args:
        prokka_file {list}: Prokka results files.
        output_file {str}: Summary output file.
    """
    
    # Get columns (because some files can have different columns)
    columns = set()
    for prokka_file in prokka_files:
        with open(prokka_file, "r") as prokka_f:
            file_columns = [line[:-1].split(":")[0] for line in prokka_f.readlines()]
            columns.update(file_columns)
    columns = ["Sample"] + list(columns)

    # Get data
    summary_df = pd.DataFrame(columns = columns)

    for prokka_file in prokka_files:
        with open(prokka_file, "r") as prokka_f:
            row = {'Sample': prokka_file.split("/")[-1].split(".")[0]}
            for line in prokka_f.readlines():
                row[line.split(":")[0]] = line[:-1].split(":")[1][1:]
        summary_df = summary_df.append(row, ignore_index = True)

    # Export to tsv
    summary_df.to_csv(output_file, sep="\t",index=False)


def refactor_gff_from_prokka(gff_input, gff_output):
    """
    Refactor gff files from prokka to an AMR-ready format.

    Arguments:
        input_gff {[type]} -- [description]
        output_gff {[type]} -- [description]
    """
    with open(gff_input) as in_file, open(gff_output, 'w+') as out_file:
        for line in in_file:
            if line.startswith("##FASTA"):
                break
            if "ID=" in line:
                line_id = line.split("\tID=")[1].split(";")[0]
                if "Name=" in line:
                    line_name = line.split("Name=")[1].split(";")[0]
                    line = line.replace("Name="+line_name, "Name="+line_id+";OldName="+line_name)
                else:
                    line = line.replace("ID="+line_id, "ID="+line_id+";Name="+line_id)
            out_file.write(line)


def dfast_call(locus_tag, contigs_file, output_dir, sample_basename, organism, threads=None):
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
                "--sort_sequence", "false",
                "--minimum_length", str(cfg.config["min_contig_len"]),
                "--use_original_name", "true",
                "--step", "1",
                "--organism", organism, 
                "--strain", sample_basename,
                "--locus_tag_prefix", locus_tag,
                "--out", output_dir]
    if threads:
        arguments.extend(["--cpu", str(threads)])

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


def refactor_gff_from_dfast(gff_input, gff_output, gff_roary_format):
    """
    Dfast sets a generic ID in its gff file, so we replace it with one related to each sample.
    
    Arguments:
        gff_input {string} -- Input gff file (from dfast)
        gff_output {string} -- Output gff file.
    """
    with open(gff_input) as in_file, open(gff_output, 'w+') as out_file, open(gff_roary_format, "w") as roary_gff:
        in_assembly = False # Flag to avoid writing the assembly part in the gffs for AMRfinder but keep it in roary ones
        for line in in_file:
            line = line.strip()
            line_w_name = ""
            if line.startswith("##FASTA"):
                in_assembly = True
            if "locus_tag=" in line:
                locus_tag = line.split("locus_tag=")[1].split(";")[0]
                line_id = line.split("\tID=")[1].split(";")[0]

                line = line.replace(line_id, locus_tag, 1)
                line_w_name = line + ";" + "Name="+line_id+"|"+locus_tag+"\n"
            if not in_assembly:
                out_file.write(line_w_name)
            roary_gff.write(line+"\n")

        

def roary_call(input_files, output_dir, campype_output_folder, threads=None):
    """
    Roary call.
    
    Arguments:
        input_files {list} -- GFF files from prokka.
        output_dir {string} -- Output directory.
    
    Returns:
        {int} -- Execution state (0 if everything is all right)
    """
    arguments = ["roary", "-f", output_dir, "-s", "-v"]
    if cfg.config["pangenome"]["no_split_paralogs"]:
        arguments.append("-s")
    if cfg.config["pangenome"]["minid"]:
        arguments.extend(["-i", str(cfg.config["pangenome"]["minid"])])
    arguments.extend(input_files)
    if threads:
        arguments.extend(["-p", str(threads)])

    ex_state = call(arguments)
    # Set Roary output directory name
    for root, dirs, _files in os.walk(campype_output_folder):
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
        general_dir {string} -- General output directory.
        input_newick {string} -- Filename (and route) to newick input file.
        input_gene_presence_absence {[type]} -- Filename (and route) to gene presence/absence input file.
        output_dir {[type]} -- Route to output files.
    
    Returns:
        {int} -- Execution state (0 if everything is all right)
    """
    arguments = ["python", "utils/roary_plots.py", input_newick, input_gene_presence_absence]
    ex_state = call(arguments)
    
    # Roary_plots saves output files in the current directory, so we move them to our own
    for filename in os.listdir("."):
        # I'm being this specific so we don't move any other similar file
        if filename.startswith("pangenome_frequency.png") or filename.startswith("pangenome_matrix.png") or filename.startswith("pangenome_pie.png"):
            shutil.move("./"+filename, output_dir+"/"+filename)
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
        call(arguments, stdout=f)
    call(["cat", bypass])

    # Create a new dataset with prinseq stats for every sample
    with open(bypass, "r") as bypass_file:
        data = csv.reader(bypass_file, delimiter="\t")

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
        call(arguments, stdout=f)
    call(["cat", bypass1])

    bypass2 = "bypass2.tsv"
    with open(bypass2, "w") as f:
        arguments = ["prinseq-lite.pl", 
                "-fastq", notcombined1, 
                "-fastq2", notcombined2, 
                "-stats_info", 
                "-stats_len"]
        call(arguments, stdout=f)
    call(["cat", bypass2])

    # Create a new dataset with prinseq stats for every sample
    with open(bypass1, "r") as bypass_file1:
        with open(bypass2, "r") as bypass_file2:
            data1 = csv.reader(bypass_file1, delimiter="\t")
            data2 = csv.reader(bypass_file2, delimiter="\t")

            if not os.path.isfile(output):
                with open(output, "w") as out_file:
                    outputwriter = csv.writer(out_file, delimiter="\t")
                    outputwriter.writerow(["Sample", "JoinReads (%)", "JoinReads", "JoinLenMeanReads", "JoinLenSdReads", "UnjoinR1Reads", "UnjoinR1LenMeanReads", "UnjoinR1LenSdReads", "UnjoinR2Reads", "UnjoinR2LenMeanReads", "UnjoinR2LenSdReads"])
                    
            
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

                data_dict["JoinReads (%)"] = str(round(data_dict["JoinReads"] / (previous_data["R1Reads"] * 2) * 100, 2))

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
                            
                outputwriter.writerow([ data_dict["Sample"], data_dict["JoinReads (%)"], data_dict["JoinReads"], data_dict["JoinLenMeanReads"], 
                                        data_dict["JoinLenSdReads"], data_dict["UnjoinR1Reads"], data_dict["UnjoinR1LenMeanReads"], data_dict["UnjoinR1LenSdReads"], 
                                        data_dict["UnjoinR2Reads"], data_dict["UnjoinR2LenMeanReads"], data_dict["UnjoinR2LenSdReads"]])
                
    call(["rm", bypass1])
    call(["rm", bypass2])
    return data_dict


def generate_report(samples, prinseq_dir, assembly_dir, annotation_dir, mauve_dir, out_dir, info_pre_QC, info_post_QC, kraken_combined_report, info_post_flash, mlst_file, vir_matrix_file, snippy_summary, custom_VFDB, amrfinder_matrix_file=False, reference_genome_basename=None):
    """
    Creates the final report.
    
    Arguments:
        samples {list} -- Samples list.
        prinseq_dir {string} -- Prinseq results directory.
        assembly_dir {string} -- Assembly results directory.
        annotation_dir {string} -- Annotator results directory.
        mauve_dir {string} --  Mauves results directory.
        out_dir {string} -- Output directory.
        infro_pre_QC {dict} -- Prinseq information before running QC.
        info_post_QC {dict} -- Prinseq information after running QC.
        kraken_combined_report {string} -- Kraken report file.
        info_post_flash {dict} -- Prinseq information after running flash.
        mlst_file {string} -- MLST results file.
        vir_matrix_file {string} -- Matrix file containing virulence genes information.
        snippy_summary {string} -- Snippy summary file.
        reference_genome_basename {string} -- Reference genome basename.
    """

    assembly_report = pd.read_csv(assembly_dir+"/"+"quality_assembly_report.tsv", sep="\t")
    if reference_genome_basename:
        ref_assembly_report = pd.read_csv(assembly_dir+"/"+reference_genome_basename+"/"+reference_genome_basename+"_assembly_statistics/"+"report.tsv", sep="\t")
    
    if cfg.config["MLST"]["run_mlst"]:
        if cfg.config["MLST"]["include_cc"]:
            mlst_data = pd.read_csv(mlst_file, sep="\t", index_col=0)
            # Sample column (index) is a string
            mlst_data.index = mlst_data.index.astype(str)
        else:
            # Read only the first 3 columns
            mlst_data = pd.read_csv(mlst_file, sep="\t", index_col=0, usecols=[0,1,2], header=None)
            mlst_data.columns = ["Genus", "ST"]
            # Rename index as Sample
            mlst_data.index.name = "Sample"
            # Sample column (index) is a string
            mlst_data.index = mlst_data.index.astype(str)
            # Filter df to keep only the basename of the index column
            mlst_data.index = mlst_data.index.str.split("/").str[-1]
            mlst_data.index = mlst_data.index.str.split(".").str[0]
            # Print index values
            print(mlst_data.index.values)


    if cfg.config["virulence_genes"]["run_virulence_genes_prediction"]:
        vir_matrix = pd.read_csv(vir_matrix_file, sep="\t", skipfooter=1, engine="python")
        vir_total_by_categories = vir_matrix.groupby(["Type"]).sum().sum(axis=1).to_dict()
        vir_types_summary = vir_matrix.groupby(["Type"]).sum().to_dict()

    fasta_mode = cfg.config["assembled_genomes"]

    if fasta_mode:    
        df_columns = [ "Sample", "ContigLen", "ST", "clonal_complex", "CDS", "CRISPRs",
                    "rRNAs", "tRNAs"]
    else:
        df_columns = [ "Sample", "Reads", "ReadLen", "ReadsQC", "ReadsQCLen", "JoinReads", "JoinReadsLen", "Contigs", 
                    "GenomeLen", "ContigLen", "N50", "GC", "DepthCov (X)", "ST", "clonal_complex", "CDS", "CRISPRs",
                    "rRNAs", "tRNAs"]
    
    if cfg.config["species_identification"]["run_species_identification"]:
        df_columns.insert(df_columns.index("ST"), "Species")
        df_columns.insert(df_columns.index("ST"), "Species_percentage")

        # Get species information (header is included in the tsv)
        kraken_report_df = pd.read_csv(kraken_combined_report, sep="\t", header=0)
        # Sample and Species columns are strings
        kraken_report_df["Sample"] = kraken_report_df["Sample"].astype(str)
        kraken_report_df["Species"] = kraken_report_df["Species"].astype(str)

    if amrfinder_matrix_file:
        amr_data = pd.read_csv(amrfinder_matrix_file, sep="\t", skipfooter=1, engine="python")
        amr_data["Sample"] = amr_data["Sample"].astype(str)
        amr_data = amr_data.set_index("Sample")
        df_columns.extend(amr_data.columns)
        

    csv_report = pd.DataFrame(columns=df_columns)

    for sample in samples+[reference_genome_basename]:
        if sample is not None:

            if sample == reference_genome_basename:
                assembly_report = ref_assembly_report
            
            # "Contigs": Number of contigs of the genome (> 500bp).
            n_contigs = int(assembly_report.loc[assembly_report['Assembly'].isin(["# contigs"])][sample])
            
            # "GenomeLen": Length (bp) of the genome.
            genome_len = int(assembly_report.loc[assembly_report['Assembly'].isin(["Total length"])][sample])

            # "N50": “Length of the smallest contig in the set that contains the fewest (largest) contigs whose combined length represents at least 50% of the assembly” (Miller et al., 2010).
            n50 = float(assembly_report.loc[assembly_report['Assembly'].isin(["N50"])][sample])
            
            # "GC": GC content (%) of the draft genome
            gc = float(assembly_report.loc[assembly_report['Assembly'].isin(["GC (%)"])][sample])

            
            if  sample != reference_genome_basename:
                if not fasta_mode:

                    # "Rads": Total number of reads after quality filtering
                    reads = info_pre_QC[sample]["R1Reads"] + info_pre_QC[sample]["R2Reads"]
                    
                    # "ReadLen": Average read length (bp) before quality control.
                    readlen = np.mean([info_pre_QC[sample]["R1LenMean"], info_pre_QC[sample]["R2LenMean"]])

                    # "ReadsQC": Total number of reads after quality control.
                    readsqc = info_post_QC[sample]["R1Reads"] + info_post_QC[sample]["R2Reads"]

                    # "ReadsQCLen": Average read length (bp) after quality control.
                    readsqclen = np.mean([info_post_QC[sample]["R1LenMean"], info_post_QC[sample]["R2LenMean"]])
                    
                    # "DepthCov (X)": Number of times each nucleotide position in the draft genome has a read that align to that position.
                    depthcov = round((info_post_QC[sample]["R1Reads"] * info_post_QC[sample]["R1LenMean"] + info_post_QC[sample]["R2Reads"] * info_post_QC[sample]["R2LenMean"])/ genome_len, 0)

                    if cfg.config["merge_reads"]:
                        # "JoinReads": Total combined reads.
                        joinreads = info_post_flash[sample]["JoinReads"]

                        # "JoinReadsLen: Mean length of combined reads.
                        joinreadslen = info_post_flash[sample]["JoinLenMeanReads"]
                        
                    else:
                        joinreads = 0
                        joinreadslen = 0

                # "ContigLen": Average contig length (bp) (> 500bp).
                contig_len_summatory = 0
                contig_counter = 0
                for record in SeqIO.parse(mauve_dir+"/"+sample+".fasta", "fasta"):
                    contig_len_summatory += len(record.seq)
                    contig_counter += 1
                avg_contig_len = contig_len_summatory/contig_counter
                
                # Kraken data
                if cfg.config["species_identification"]["run_species_identification"]:
                    # "Species": Species name.
                    species = kraken_report_df.loc[kraken_report_df["Sample"] == sample]["Species"].values[0]
                    
                    # "Species_percentage": Percentage of species.
                    species_percentage = kraken_report_df.loc[kraken_report_df["Sample"] == sample]["Species_percentage"].values[0]
            

            else:
                # "ContigLen": Average contig length (bp) (> 500bp).
                contig_len_summatory = 0
                contig_counter = 0
                for record in SeqIO.parse(cfg.config["reference_genome"]["file"], "fasta"):
                    contig_len_summatory += len(record.seq)
                    contig_counter += 1
                avg_contig_len = contig_len_summatory/contig_counter
                if not fasta_mode:
                    reads = readlen = readsqc = readsqclen = joinreads = joinreadslen = depthcov = 0
                
                # Kraken data for reference genome
                if cfg.config["species_identification"]["run_species_identification"]:
                    species = ""
                    species_percentage = ""

            # Column ST (MLST)
            if cfg.config["MLST"]["run_mlst"]:
                if cfg.config["MLST"]["include_cc"]:
                    st = mlst_data.loc[sample]["ST"]
                    # Column clonal_complex (MLST)
                    clonal_complex = mlst_data.loc[sample]["clonal_complex"]
                else:
                    st = mlst_data.loc[sample]["ST"]
                    clonal_complex = None
            else:
                st = None
                clonal_complex = None

            
            
            cds = "0"
            crisprs = "0"
            rrnas = "0"
            trnas = "0"
            hy_prot = 0
            
            if cfg.config["annotation"]["run_annotation"]:
                if cfg.config["annotation"]["annotator"] == "prokka":
                    
                    with open(annotation_dir+"/"+sample+"/"+sample+".txt") as stats_file:
                        for line in stats_file:
                            # Column 'CDS'. 
                            if line.startswith("CDS:"):
                                cds = line.split()[-1].replace("\n", "")

                            # Column 'CRISPRs'. 
                            elif line.startswith("repeat_region:") or line.startswith("CRISPR:"):
                                crisprs = line.split()[-1].replace("\n", "")

                            # Column 'rRNAs'. 
                            elif line.startswith("rRNA:"):
                                rrnas = line.split()[-1].replace("\n", "")

                            # Column 'tRNAs'. 
                            elif line.startswith("tRNA:"):
                                trnas = line.split()[-1].replace("\n", "")

                    with open(annotation_dir+"/"+sample+"/"+sample+".faa") as faa_file:
                        # Column 'Hypothetical_proteins'
                        hy_prot = faa_file.read().count("hypothetical protein")


                elif cfg.config["annotation"]["annotator"] == "dfast":
                    with open(annotation_dir+"/"+sample+"/"+sample+"_statistics.txt") as stats_file:
                        for line in stats_file:
                            # Column 'CDS'.
                            if line.startswith("Number of CDSs"):
                                cds = line.split("\t")[-1].replace("\n", "")
                            
                            # Column 'CRISPRs'.
                            elif line.startswith("Number of CRISPRs"):
                                crisprs = line.split("\t")[-1].replace("\n", "")

                            # Column 'rRNAs'. 
                            elif line.startswith("Number of rRNAs"):    
                                rrnas = line.split("\t")[-1].replace("\n", "")

                            # Column 'tRNAs'. 
                            elif line.startswith("Number of tRNAs"):
                                trnas = line.split("\t")[-1].replace("\n", "")
                    # TODO Column 'Hypothetical proteins'
                    with open(annotation_dir+"/"+sample+"/"+sample+".faa") as faa_file:
                        # Column 'Hypothetical_proteins'
                        hy_prot = faa_file.read().count("hypothetical protein")

            if fasta_mode:  
                report_dict = { "Sample": sample, 
                            "Contigs": n_contigs, 
                            "GenomeLen": genome_len, 
                            "ContigLen": round(avg_contig_len, 2), 
                            "N50": round(n50, 0),
                            "GC": round(gc, 2),
                            "ST": st,
                            "clonal_complex": clonal_complex,
                            "CDS": cds,
                            "Hypothetical proteins": round(hy_prot, 0),
                            "CRISPRs": crisprs,
                            "rRNAs": rrnas,
                            "tRNAs": trnas}
            else:
                report_dict = { "Sample": sample, 
                            "Reads": round(reads, 0), 
                            "ReadLen": readlen,
                            "ReadsQC": readsqc,
                            "ReadsQCLen": readsqclen, 
                            "JoinReads": round(joinreads, 2), 
                            "JoinReadsLen": joinreadslen, 
                            "Contigs": n_contigs, 
                            "GenomeLen": genome_len, 
                            "ContigLen": round(avg_contig_len, 2), 
                            "N50": round(n50, 0),
                            "GC": round(gc, 2),
                            "DepthCov (X)": round(depthcov, 2),
                            "ST": st,
                            "clonal_complex": clonal_complex,
                            "CDS": cds,
                            "Hypothetical proteins": round(hy_prot, 0),
                            "CRISPRs": crisprs,
                            "rRNAs": rrnas,
                            "tRNAs": trnas}
            
            # Kraken
            if cfg.config["species_identification"]["run_species_identification"]:
                report_dict["Species"] = species
                report_dict["Species_percentage"] = species_percentage

            # Columna 20 y siguientes. Se obtendrán a partir del archivo AMR_genes_AMRFinder_matrix.txt. Se copiarán las columnas 2 y siguientes, de forma que se corresponda la información de cada cepa.
            if amrfinder_matrix_file:
                report_dict.update(amr_data.loc[sample].to_dict())
            
            # Column blocks to sort columns in final file
            df_1st_column_block = [*report_dict.keys()]
            df_2nd_column_block = []
            
            occurrences = 0     # Total number of occurrences of that category in Campylobacter_custom_VFDB.txt
            total_occurrences = 0
            if cfg.config["virulence_genes"]["run_virulence_genes_prediction"]:
                # Information from VF_matrix
                # For each sample its column has the sum of genes present from each category
                virulence_dict = {}
                for category, _ in vir_total_by_categories.items():
                    with open(custom_VFDB) as custom_DB:
                        occurrences = sum(category.lower() in line.split(" ")[1].lower() for line in custom_DB if line.startswith(">"))
                        total_occurrences += occurrences
                    virulence_dict[category+"("+str(occurrences)+")"] = vir_types_summary[sample][category]
                df_2nd_column_block = [*virulence_dict.keys()]     # Column blocks to sort columns in final file
                virulence_dict["Total ("+str(total_occurrences)+")"] = vir_matrix.groupby(["Type"]).sum()[sample].sum()
                
                report_dict.update(virulence_dict)


            if cfg.config["run_variant_calling"] and cfg.config["reference_genome"]["file"]:
                #Snippy block 
                if sample != reference_genome_basename:
                    snippy_data = pd.read_csv(snippy_summary, sep="\t", dtype={'Sample':str})
                    snippy_data.set_index("Sample", inplace=True)
                    df_3rd_column_block = [*snippy_data.columns]
                    report_dict.update(snippy_data.loc[sample].to_dict())
            else:
                df_3rd_column_block = None


            csv_report = csv_report.append(report_dict, ignore_index=True)
        
    new_colums_order = [*df_1st_column_block, *df_2nd_column_block, "Total ("+str(total_occurrences)+")"]
    if df_3rd_column_block:
        new_colums_order += df_3rd_column_block
    
    csv_report = csv_report.reindex(columns=new_colums_order)
    csv_report.to_csv(out_dir+"/campype_report.csv", sep="\t", index=False)


if __name__ == "__main__":

    # Welcome
    welcome("resources/campype_ascii.txt")

    # Get config file parameters
    annotator = cfg.config["annotation"]["annotator"]
    fasta_mode = cfg.config["assembled_genomes"]
    # Create output directories
    now = datetime.datetime.now()

    # Get reference files from campype_config.py
    adapters_file =  cfg.config["trim_adaptors"]["adapters_reference_file"]
    reference_genome_file = cfg.config["reference_genome"]["file"]
    proteins_file = cfg.config["virulence_genes"]["blast"]["proteins_reference_file"]
    
    output_folder = sys.argv[1]

    # Generate json config file from campype_config.py
    with open(output_folder+"/campype_config.json", 'w') as json_file:
        json.dump(cfg.config, json_file)

    trimmomatic_dir = output_folder+"/tmp_Trimmomatic_filtering"
    if cfg.config["trim_adaptors"]["run_trim_adaptors"]:
        prinseq_dir = output_folder+"/Trimmomatic_and_Prinseq_filtering"
    else:
        prinseq_dir = output_folder+"/Prinseq_filtering"
    kraken_dir = output_folder+"/Kraken_identification"
    kraken_unified_report = kraken_dir+"/bacteria_identification.tab"
    flash_dir = output_folder+"/Flash_read_extension"
    spades_dir = output_folder+"/SPAdes_assembly"
    contigs_dir = output_folder+"/Draft_genomes"
    mauve_dir = output_folder+"/Mauve_reordered_draft_genomes"
    snps_dir = output_folder+"/SNP_SNIPPY"
    mlst_dir = output_folder+"/MLST"
    vir_dir = output_folder+"/Virulence_genes"
    plasmid_dir = output_folder+"/Plasmids"
    amr_analysis_dir = output_folder+"/Antimicrobial_resistance_genes"
    amr_analysis_dir_abr = amr_analysis_dir+"/ABRicate"
    amr_analysis_dir_amrfinder = amr_analysis_dir+"/AMRFinder"
    prokka_dir = output_folder+"/Prokka_annotation"

    annotation_dir = "" # Will be defined later
    proteins_database_name = ""
    
    dfast_dir = output_folder+"/Dfast_annotation"
    roary_dir = output_folder+"/Roary_pangenome"
    roary_input_dir = roary_dir+"/input_gffs"
    roary_plots_dir = roary_dir+"/Roary_plots"
    dfast_refactor_dir = dfast_dir+"/refactored_gff"
    prokka_refactor_dir = prokka_dir+"/AMR_format_files"
    blast_proteins_dir = vir_dir+"/BLAST_inhouse_virulence_genes"
    dna_database_blast = blast_proteins_dir+"/DNA_database"

    # Create directories

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    
    if cfg.config["trim_adaptors"]["run_trim_adaptors"]:
        os.mkdir(trimmomatic_dir)

    if cfg.config["plasmids"]["run_plasmid_prediction"]:
        os.mkdir(plasmid_dir)

    os.mkdir(prinseq_dir)
    if cfg.config["species_identification"]["run_species_identification"]:
        os.mkdir(kraken_dir)
    os.mkdir(spades_dir)
    os.mkdir(contigs_dir)
    os.mkdir(mauve_dir)
    os.mkdir(snps_dir)
    os.mkdir(mlst_dir)
    
    if cfg.config["merge_reads"]:
        os.mkdir(flash_dir)

    if cfg.config["virulence_genes"]["run_virulence_genes_prediction"]:
        os.mkdir(vir_dir)
        if "blast" in cfg.config["virulence_genes"]["virulence_genes_predictor_tool"]:
            os.mkdir(blast_proteins_dir)
        
    
    if cfg.config["antimicrobial_resistance_genes"]["run_antimicrobial_resistance_genes_prediction"]:
        os.mkdir(amr_analysis_dir)
        if "abricate" in cfg.config["antimicrobial_resistance_genes"]["antimicrobial_resistance_genes_predictor_tool"]:
            os.mkdir(amr_analysis_dir_abr)
        if "amrfinder" in cfg.config["antimicrobial_resistance_genes"]["antimicrobial_resistance_genes_predictor_tool"]:
            os.mkdir(amr_analysis_dir_amrfinder)

    os.mkdir(roary_dir)
    if cfg.config["annotation"]["run_annotation"]:
        if annotator == "dfast":
            os.mkdir(dfast_dir)
            os.mkdir(dfast_refactor_dir)
            os.mkdir(roary_input_dir)
        elif annotator == "prokka":
            os.mkdir(prokka_dir)
            if cfg.config["antimicrobial_resistance_genes"]["run_antimicrobial_resistance_genes_prediction"]:
                if "amrfinder" in cfg.config["antimicrobial_resistance_genes"]["antimicrobial_resistance_genes_predictor_tool"]:
                    os.mkdir(prokka_refactor_dir)

    roary_input_files = []
    snippy_output_files = []
    summary_pre_qc = {}
    summary_post_qc = {}
    summary_post_flash = {}

    compressed_mode = False
    decompressed_samples_fw = dict()
    decompressed_samples_rv = dict()
    decompressed_samples = dict()

    # Number of threads
    n_threads = cfg.config["n_threads"]

    # Check if any of the input files does not have species/genus information
    abort_flag = False
    for sample_basename, data in read_input_files("input_files.csv"):
        if not data["Genus"] or not data["Species"]:
            if not cfg.config["species_identification"]["run_species_identification"]:
                print("ERROR: Please, fill in the genus and species information for sample", sample_basename, "or set species_identification parameter to True in campype_config.py")
                abort_flag = True
    if abort_flag:
        sys.exit(1)

    # Checking correct input files format in fasta mode
    for sample_basename, data in read_input_files("input_files.csv"):
        sample_fw = data["FW"]
        sample_rv = data["RV"]
        
        if fasta_mode:
            if sample_fw.split(".")[-1] == "gz":
                if not sample_fw.split(".")[-2] in ["fasta", "fa", "fna"]:
                    print("WRONG INPUT FORMAT:", sample_fw, flush=True)
                    print("Input must be a valid fasta format file.", flush=True)
                    sys.exit(1)
                if sample_fw.split(".")[-2] in ["fasta", "fa", "fna"]:
                    # Mauve needs fasta extension to work, we will copy this files to a new folder with the correct extension
                    call("cp " + sample_fw + " " + output_folder + "/" + sample_fw.split("/")[-1].split(".")[0] + ".fasta", shell=True)
                    compressed_mode = True # Not exactly compressed mode but we will use this mode to remove the files at the end of the program
                                        
                # Decompress files
                compressed_mode = True
                if sample_fw.split(".")[-2] == "tar":
                    # Unzip from .tar.gz
                    call("tar -xzf " + sample_fw + " -C " + output_folder, shell=True)
                    call("tar -xzf " + sample_rv + " -C " + output_folder, shell=True)
                    # Get file basename from path
                    decompressed_samples[sample_basename] = output_folder + "/" + sample_fw.split("/")[-1].split(".")[0] + ".fasta"
                else:
                    # Unzip from .gz
                    with open(output_folder + "/" + sample_fw.split("/")[-1].split(".")[0] + ".fasta", "w") as fw_file:
                        call("gunzip -c " + sample_fw, stdout=fw_file, shell=True)
                    decompressed_samples[sample_basename] = output_folder + "/" + sample_fw.split("/")[-1].split(".")[0] + ".fasta"
            else:
                if not sample_fw.split(".")[-1] in ["fasta", "fa", "fna"]:
                    print("WRONG INPUT FORMAT:", sample_fw, flush=True)
                    print("Input must be a valid fasta format file.", flush=True)
                    sys.exit(1)
        else:
            if sample_fw.split(".")[-1] == "gz":
                compressed_mode = True
                if sample_fw.split(".")[-2] == "tar":
                    # Unzip from .tar.gz
                    call("tar -xzf " + sample_fw + " -C " + output_folder, shell=True)
                    call("tar -xzf " + sample_rv + " -C " + output_folder, shell=True)
                    # Get file basename from path
                    decompressed_samples_fw[sample_basename] = output_folder + "/" + sample_fw.split("/")[-1].split(".")[0] + ".fastq"
                    decompressed_samples_rv[sample_basename] = output_folder + "/" + sample_rv.split("/")[-1].split(".")[0] + ".fastq"
                else:
                    # Unzip from .gz
                    with open(output_folder + "/" + sample_fw.split("/")[-1].split(".")[0] + ".fastq", "w") as fw_file:
                        call("gunzip -c " + sample_fw, stdout=fw_file, shell=True)
                    with open(output_folder + "/" + sample_rv.split("/")[-1].split(".")[0] + ".fastq", "w") as rv_file:
                        call("gunzip -c " + sample_rv, stdout=rv_file, shell=True)
                    decompressed_samples_fw[sample_basename] = output_folder + "/" + sample_fw.split("/")[-1].split(".")[0] + ".fastq"
                    decompressed_samples_rv[sample_basename] = output_folder + "/" + sample_rv.split("/")[-1].split(".")[0] + ".fastq"
                

    # Reference file processing
    if reference_genome_file:
        reference_genome_filename = reference_genome_file.split("/")[-1]
        reference_genome_basename = reference_genome_filename.split(".")[-2]
        
        print(Banner(f"\nProcessing reference sequence ({reference_genome_basename}) before your genomes of interest  \n"), flush=True)

        # Quast call for reference file

        # Create Quast output directories
        quast_sample_dir = reference_genome_basename+"_assembly_statistics"
        if fasta_mode: # In fasta mode the Spades directory does not exist
            quast_dir = output_folder+"/Quast/"
            assembly_dir = quast_dir
            os.makedirs(assembly_dir+quast_sample_dir)
        else:
            assembly_dir = spades_dir
            quast_dir = assembly_dir+"/"+reference_genome_basename+"/"
            os.mkdir(quast_dir)
            os.mkdir(quast_dir+quast_sample_dir)


        print(Banner(f"\nQuast on reference sequence ({reference_genome_basename})\n"), flush=True)
        quast_call( input_file=reference_genome_file,
                    output_dir=assembly_dir+"/"+reference_genome_basename+"/"+quast_sample_dir,
                    min_contig_len=cfg.config["min_contig_len"],
                    threads = n_threads)

        if cfg.config["annotation"]["run_annotation"]:    
            # Annotate reference fasta file 
            if annotator.lower() == "dfast":
                # Dfast call
                annotation_dir = dfast_dir
                print(Banner(f"\nAnnotating reference sequence: Dfast\n"), flush=True)
                dfast_call( locus_tag=reference_genome_file+"_L",
                            contigs_file=reference_genome_file,
                            output_dir=dfast_dir+"/"+reference_genome_basename,
                            sample_basename=reference_genome_basename,
                            organism=cfg.config["reference_genome"]["genus"]+" "+cfg.config["reference_genome"]["species"],
                            threads = n_threads)
                refactor_gff_from_dfast(dfast_dir+"/"+reference_genome_basename+"/"+reference_genome_basename+".gff",
                                        dfast_refactor_dir+"/+"+reference_genome_basename+".gff",
                                        roary_input_dir+"/+"+reference_genome_basename+".gff")
                # Set roary input files (renaming to get reference file first)
                # os.rename(dfast_refactor_dir+"/"+reference_genome_basename+".gff", dfast_refactor_dir+"/+"+reference_genome_basename+".gff")
                roary_input_files.append(roary_input_dir+"/+"+reference_genome_basename+".gff")
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
                            rawproduct=cfg.config["annotation"]["prokka"]["rawproduct"],
                            threads = n_threads
                            )
                # Set roary input files (renaming to get reference file first)
                os.rename(annotation_dir+"/"+reference_genome_basename+"/"+reference_genome_basename+".gff",
                        annotation_dir+"/"+reference_genome_basename+"/+"+reference_genome_basename+".gff")
                if cfg.config["antimicrobial_resistance_genes"]["run_antimicrobial_resistance_genes_prediction"]:
                    if "amrfinder" in cfg.config["antimicrobial_resistance_genes"]["antimicrobial_resistance_genes_predictor_tool"]:
                        refactor_gff_from_prokka(prokka_dir+"/"+reference_genome_basename+"/+"+reference_genome_basename+".gff", prokka_refactor_dir+"/"+"+"+reference_genome_basename+".gff")
                roary_input_files.append(annotation_dir+"/"+reference_genome_basename+"/+"+reference_genome_basename+".gff")

    
    # Workflow Starts (for standard samples)
    sample_counter = 0
    samples_basenames = [] # Keeping track of them
    kraken_reports = []
    n_samples = len(read_input_files("input_files.csv"))

    for sample_basename, data in read_input_files("input_files.csv"):
        
        sample_fw = data["FW"]
        sample_rv = data["RV"]
        genus = data["Genus"]
        species = data["Species"]

        step_counter = 1 # Just to let the user know the number of each step
        sample_counter += 1
        samples_basenames.append(sample_basename)
        
        if not fasta_mode:
            # Run trimmomatic or not
            if cfg.config["trim_adaptors"]["run_trim_adaptors"]:
                
                print(Banner(f"\nStep {step_counter} for sequence {sample_counter}/{n_samples} ({sample_basename}): Trimmomatic\n"), flush=True)
                trimmomatic_call(input_file1=sample_fw,
                                input_file2=sample_rv,
                                phred="-phred33",
                                trimfile="ILLUMINACLIP:"+adapters_file+":1:30:11",
                                paired_out_file1=trimmomatic_dir+"/"+sample_basename+"_R1_paired.fastq",
                                unpaired_out_file1=trimmomatic_dir+"/"+sample_basename+"_R1_unpaired.fastq",
                                paired_out_file2=trimmomatic_dir+"/"+sample_basename+"_R2_paired.fastq",
                                unpaired_out_file2=trimmomatic_dir+"/"+sample_basename+"_R2_unpaired.fastq",
                                threads = n_threads)
                step_counter += 1
                prinseq_input1 = trimmomatic_dir+"/"+sample_basename+"_R1_paired.fastq"
                prinseq_input2 = trimmomatic_dir+"/"+sample_basename+"_R2_paired.fastq"
            else:
                prinseq_input1 = decompressed_samples_fw[sample_basename]
                prinseq_input2 = decompressed_samples_rv[sample_basename]


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
            print("\nRefactoring prinseq output files\n", flush=True)
            prinseq_files = refactor_prinseq_output(prinseq_dir+"/"+sample_basename, sample_basename)

            # Kraken call
            if cfg.config["species_identification"]["run_species_identification"]:
                print(Banner(f"\nStep {step_counter} for sequence {sample_counter}/{n_samples} ({sample_basename}): Kraken\n"), flush=True)
                
                kraken_db_file = cfg.config["species_identification"]["species_identification_database"]
                
                # Check if database exists
                if not os.path.exists(kraken_db_file):
                    print("ERROR: Kraken database file does not exist.", flush=True)
                    sys.exit(1)
                kraken_out_report_file = kraken_dir+"/"+sample_basename+".tab"
                kraken_genus, kraken_species = kraken_call(db_file=kraken_db_file, 
                                                           output_file=kraken_out_report_file, 
                                                           fw_file=prinseq_files["R1"], 
                                                           rv_file=prinseq_files["R2"],
                                                           threads = n_threads)
                kraken_reports.append((sample_basename, kraken_out_report_file))
                if not genus or not species:
                    genus = kraken_genus
                    species = kraken_species
                step_counter += 1


            # Flash call
            if cfg.config["merge_reads"]:
                print(Banner(f"\nStep {step_counter} for sequence {sample_counter}/{n_samples} ({sample_basename}): Flash\n"), flush=True)
                flash_call(input_file_1=prinseq_files["R1"],
                        input_file_2=prinseq_files["R2"],
                        output_filename=sample_basename,
                        output_dir=flash_dir+"/"+sample_basename,
                        threads = n_threads)
                step_counter += 1

            # Quality reports
            print(Banner(f"\nStep {step_counter} for sequence {sample_counter}/{n_samples} ({sample_basename}): Read statistics\n"), flush=True)
            step_counter += 1
            
            report_pre_qc = get_reads_table(decompressed_samples_fw[sample_basename], decompressed_samples_rv[sample_basename], sample_basename, prinseq_dir+"/reads_statistics_beforeQC.tsv", False)

            report_post_qc = get_reads_table(prinseq_files["R1"], prinseq_files["R2"], sample_basename, prinseq_dir+"/reads_statistics_afterQC.tsv", True)
            if cfg.config["merge_reads"]:
                report_post_flash = get_flash_reads_table(flash_dir+"/"+sample_basename+"/"+sample_basename+".extendedFrags.fastq", 
                                    flash_dir+"/"+sample_basename+"/"+sample_basename+".notCombined_1.fastq",
                                    flash_dir+"/"+sample_basename+"/"+sample_basename+".notCombined_2.fastq",
                                    sample_basename, flash_dir+"/reads_statistics_FLASH.tsv", report_post_qc)

            summary_pre_qc[sample_basename] = report_pre_qc
            summary_post_qc[sample_basename] = report_post_qc
            if cfg.config["merge_reads"]:
                summary_post_flash[sample_basename] = report_post_flash

            # Create SPAdes output directories
            os.mkdir(spades_dir+"/"+sample_basename)

            # SPAdes call
            print(Banner(f"\nStep {step_counter} for sequence {sample_counter}/{n_samples} ({sample_basename}): SPAdes\n"), flush=True)
            if cfg.config["merge_reads"]:        
                spades_call(forward_sample=flash_dir+"/"+sample_basename+"/"+sample_basename+".notCombined_1.fastq",
                            reverse_sample=flash_dir+"/"+sample_basename+"/"+sample_basename+".notCombined_2.fastq",
                            sample=sample_basename,
                            out_dir=spades_dir,
                            merged_sample=flash_dir+"/"+sample_basename+"/"+sample_basename+".extendedFrags.fastq",
                            threads = n_threads)
            else:
                spades_call(forward_sample=prinseq_dir+"/"+sample_basename+"/"+sample_basename+"_R1.fastq",
                            reverse_sample=prinseq_dir+"/"+sample_basename+"/"+sample_basename+"_R2.fastq",
                            sample=sample_basename,
                            out_dir=spades_dir,
                            threads = n_threads)
            step_counter += 1

            # Get minimum contig length
            # min_contig_threshold = get_reads_length(prinseq_dir+"/"+sample_basename+"/"+sample_basename+"_R1.fastq") * 2
            min_contig_threshold = cfg.config["min_contig_len"]

            # Trim short contigs and shorten sequences id
            contigs_trim_and_rename(contigs_file=spades_dir+"/"+sample_basename+"/"+"contigs.fasta",
                                    output_filename=sample_basename+".fasta",
                                    output_dir=contigs_dir,
                                    min_len=min_contig_threshold)
            sample_contigs = contigs_dir+"/"+sample_basename+".fasta"
        else:
            min_contig_threshold = cfg.config["min_contig_len"]
            # If in fasta mode, we just skip everything until this point
            if compressed_mode:
                sample_contigs = decompressed_samples[sample_basename]
            else:
                sample_contigs = data["FW"]

            # Kraken call
            if cfg.config["species_identification"]["run_species_identification"]:
                print(Banner(f"\nStep {step_counter} for sequence {sample_counter}/{n_samples} ({sample_basename}): Kraken\n"), flush=True)
                
                kraken_db_file = "./db/minikraken-DB/minikraken_8GB_20200312"
                
                # Check if database exists
                if not os.path.exists(kraken_db_file):
                    print("ERROR: Kraken database file does not exist.", flush=True)
                    sys.exit(1)

                kraken_out_report_file = kraken_dir+"/"+sample_basename+".tab" 
                kraken_genus, kraken_species = kraken_call(db_file=kraken_db_file, 
                                                           output_file=kraken_out_report_file, 
                                                           fw_file=sample_fw,
                                                           threads = n_threads)
                if not genus or not species:
                    genus = kraken_genus
                    species = kraken_species
                kraken_reports.append((sample_basename, kraken_out_report_file))
                step_counter += 1
            
        # Reordering contigs by a reference genome with MauveCM
        if cfg.config["reference_genome"]["file"]:
            print(Banner(f"\nStep {step_counter} for sequence {sample_counter}/{n_samples} ({sample_basename}): reordering {sample_basename} genome against reference genome\n"), flush=True)
            draft_contigs = mauve_call(output_folder=mauve_dir,
                                        reference_sequence=reference_genome_file,
                                        input_contigs=sample_contigs,
                                        sample_basename=sample_basename)
            step_counter += 1
            draft_contigs_dir = mauve_dir
        else:
            if fasta_mode:
                draft_contigs = sample_contigs
                # Copy contigs to contigs_dir
                shutil.copyfile(draft_contigs, contigs_dir+"/"+sample_basename+".fasta")
            else:
                draft_contigs = contigs_dir+"/"+sample_basename+".fasta"
            draft_contigs_dir = contigs_dir


        # Create Quast output directories
        quast_sample_dir = sample_basename+"_assembly_statistics"
        if fasta_mode: # In fasta mode the Spades directory does not exist
            quast_dir = output_folder+"/Quast/"
            assembly_dir = quast_dir
            os.makedirs(assembly_dir+quast_sample_dir)
        else:
            assembly_dir = spades_dir
            quast_dir = assembly_dir+"/"+sample_basename+"/"
            os.mkdir(quast_dir+quast_sample_dir)
            
            
        

        # Quast call
        print(Banner(f"\nStep {step_counter} for sequence {sample_counter}/{n_samples} ({sample_basename}): Quast\n"), flush=True)
        quast_call( input_file=draft_contigs,
                    output_dir=assembly_dir+"/"+sample_basename+"/"+quast_sample_dir,
                    min_contig_len=min_contig_threshold,
                    threads = n_threads)
        step_counter += 1

        if cfg.config["annotation"]["run_annotation"]:
            # Annotation (Prokka or dfast)
            if annotator.lower() == "dfast":
                # Dfast call
                annotation_dir = dfast_dir
                print(Banner(f"\nStep {step_counter} for sequence {sample_counter}/{n_samples} ({sample_basename}): Dfast\n"), flush=True)
                dfast_call(locus_tag=sample_basename+"_L",
                        contigs_file=draft_contigs,
                        output_dir=dfast_dir+"/"+sample_basename,
                        sample_basename=sample_basename,
                        organism=genus+" "+species,
                    threads = n_threads)
                step_counter += 1

                refactor_gff_from_dfast(dfast_dir+"/"+sample_basename+"/"+sample_basename+".gff", 
                                        dfast_refactor_dir+"/"+sample_basename+".gff",
                                        roary_input_dir+"/"+sample_basename+".gff")
                # Set roary input files
                roary_input_files.append(roary_input_dir+"/"+sample_basename+".gff")

            else:
                # Prokka call
                annotation_dir = prokka_dir
                print(Banner(f"\nStep {step_counter} for sequence {sample_counter}/{n_samples} ({sample_basename}): Prokka\n"), flush=True)
                prokka_call(locus_tag=sample_basename+"_L",
                            output_dir=prokka_dir+"/"+sample_basename,
                            prefix=sample_basename,
                            input_file=draft_contigs,
                            genus=genus,
                            species=species,
                            strain=sample_basename,
                            proteins=cfg.config["reference_genome"]["proteins"],
                            rawproduct=cfg.config["annotation"]["prokka"]["rawproduct"],
                            threads = n_threads)
                step_counter += 1
                if cfg.config["antimicrobial_resistance_genes"]["run_antimicrobial_resistance_genes_prediction"]:
                    if "amrfinder" in cfg.config["antimicrobial_resistance_genes"]["antimicrobial_resistance_genes_predictor_tool"]:
                        refactor_gff_from_prokka(prokka_dir+"/"+sample_basename+"/"+sample_basename+".gff", prokka_refactor_dir+"/"+sample_basename+".gff")
                # Set roary input files
                roary_input_files.append(annotation_dir+"/"+sample_basename+"/"+sample_basename+".gff")
        

        if cfg.config["run_variant_calling"]:
            print(Banner(f"\nStep {step_counter} for sequence {sample_counter}/{n_samples} ({sample_basename}): SNIPPY\n"), flush=True)
            if reference_genome_file:
                # SNPs identification (SNIPPY)
                step_counter += 1
                reference_genome_filename = reference_genome_file.split("/")[-1]
                reference_genome_basename = reference_genome_filename.split(".")[-2]
                if cfg.config["reference_genome"]["proteins"]:
                    ref_genome = cfg.config["reference_genome"]["proteins"]
                else:
                    ref_genome = cfg.config["reference_genome"]["file"]
                snippy_call(reference_genome=ref_genome,
                            contigs=draft_contigs,
                            output_dir=snps_dir+"/"+sample_basename,
                            prefix=sample_basename,
                            threads = n_threads)
                snippy_output_files.append(snps_dir+"/"+sample_basename+"/"+sample_basename+".txt")
            else:
                # We cannot run SNIPPY without a reference genome
                print("ERROR: You must provide a reference genome to run SNIPPY.", flush=True)

    # Kraken combined report
    if cfg.config["species_identification"]["run_species_identification"]:
        kraken_report_unification(kraken_reports, kraken_unified_report)

    # Prokka summary
    if cfg.config["annotation"]["run_annotation"] and cfg.config["annotation"]["annotator"] == "prokka":
        prokka_summary_outfile = prokka_dir+"/Prokka_summary.tsv"
        prokka_summary_input_files = [prokka_dir+"/"+sample+"/"+sample+".txt" for sample in samples_basenames]
        prokka_summary(prokka_summary_input_files, prokka_summary_outfile)

    # Snippy summary
    if cfg.config["run_variant_calling"] and reference_genome_file:
        snippy_summary_outfile = snps_dir+"/Genomic_variants_summary.tsv"
        snippy_summary(snippy_output_files, snippy_summary_outfile)
    else:
        snippy_summary_outfile = None

    # Quast report unification
    quast_report_unification(assembly_dir, samples_basenames, assembly_dir)

    # MLST call
    if cfg.config["MLST"]["run_mlst"]:
        mlst_out_file = "MLST.txt"
        print(Banner(f"\nStep {step_counter}: MLST\n"), flush=True)
        step_counter += 1
        mlst_call(input_dir=draft_contigs_dir,  
                reference_file=reference_genome_file,
                output_dir=mlst_dir,
                output_filename=mlst_out_file,
                threads = n_threads)
        if cfg.config["MLST"]["include_cc"]:
            # MLST postprocessing
            mlst_postprocessing(mlst_dir+"/"+mlst_out_file, mlst_dir+"/MLST_and_CC.txt")
            mlst_data_file = mlst_dir+"/MLST_and_CC.txt"
        else:
            mlst_data_file = mlst_dir+"/"+mlst_out_file

    # Virulence genes
    if cfg.config["virulence_genes"]["run_virulence_genes_prediction"]:

        # ABRicate call (virulence genes)
        if "abricate" in cfg.config["virulence_genes"]["virulence_genes_predictor_tool"]:
            global_vf_output_file = "Virulence_genes_ABRicate.tsv"
            global_vf_matrix_file = "Virulence_genes_ABRicate_matrix.tsv"
            for vf_database in cfg.config["virulence_genes"]["abricate"]["virulence_factors_databases"]:
                print(Banner(f"\nStep {step_counter}: Virulence genes (ABRicate: "+vf_database+")\n"), flush=True)
                step_counter += 1
                vf_output_file = "Virulence_genes_ABRicate_"+vf_database+".tsv"
                vf_matrix_file = "Virulence_genes_ABRicate_"+vf_database+"_matrix.tsv"
                abricate_samples_basenames = samples_basenames.copy()
                if cfg.config["reference_genome"]["strain"]:
                    abricate_samples_basenames.append(cfg.config["reference_genome"]["strain"])        
                abricate_call(input_dir=draft_contigs_dir,
                            output_dir=vir_dir,
                            output_filename=vf_output_file,
                            database=vf_database,
                            mincov=cfg.config["virulence_genes"]["abricate"]["mincov"],
                            minid=cfg.config["virulence_genes"]["abricate"]["minid"],
                            gene_matrix_file=vf_matrix_file,
                            samples=abricate_samples_basenames,
                            threads = n_threads)
                
                # Concatenate every ABRicate output in a single file
                with open(vir_dir+"/"+global_vf_output_file, "a") as global_file, open(vir_dir+"/"+vf_output_file, "r") as current_file:
                    if os.stat(vir_dir+"/"+global_vf_output_file).st_size == 0:   # If global file is empty, i.e. there is no header
                        global_file.write(current_file.read())
                    else:
                        for line in current_file.readlines()[1:]:
                            global_file.write(line)
                    
                with open(vir_dir+"/"+global_vf_matrix_file, "a") as global_matrix, open(vir_dir+"/"+vf_matrix_file, "r") as current_matrix:
                    if os.stat(vir_dir+"/"+global_vf_matrix_file).st_size == 0:   # If global file is empty, i.e. there is no header
                        global_matrix.write(current_matrix.read())
                    else:
                        for line in current_matrix.readlines()[1:]:
                            global_matrix.write(line)

                # Remove innecesary single files
                os.remove(vir_dir+"/"+vf_output_file)
                os.remove(vir_dir+"/"+vf_matrix_file)
                
                # End line
                with open(vir_dir+"/"+global_vf_matrix_file, "a") as matrix_file:
                    matrix_file.write("Coverage >= " + 
                                      str(cfg.config["virulence_genes"]["abricate"]["mincov"]) +
                                      " % and identity >= " +
                                      str(cfg.config["virulence_genes"]["abricate"]["minid"]) + 
                                      " % on each sample for considering a virulence gene as present.")

        # Blast call (Virulence genes)
        proteins_database_name = "inhouse_VFDB.txt"    # This is an output file name
        if "blast" in cfg.config["virulence_genes"]["virulence_genes_predictor_tool"]:
            print(Banner(f"\nStep {step_counter}: Virulence genes (BLAST against inhouse database)\n"), flush=True)
            step_counter += 1
            # Every fasta file in draft_contigs_dir is a sample
            contig_files = ([os.path.join(draft_contigs_dir, f) for f in os.listdir(draft_contigs_dir) if f.endswith(".fasta")])

            blast_output_name = "BLAST_inhouse_VFDB.txt"
            proteins_file = cfg.config["virulence_genes"]["blast"]["proteins_reference_file"]
            dna_database_blast = blast_proteins_dir+"/DNA_database"
            
            if not os.path.exists(blast_proteins_dir):
                os.mkdir(blast_proteins_dir)
            if not os.path.exists(dna_database_blast):
                os.mkdir(dna_database_blast)

            blast_contigs_files_paths = contig_files.copy()
            blast_samples_basenames = samples_basenames.copy()
            if cfg.config["reference_genome"]["file"]:
                ref_genome = cfg.config["reference_genome"]["file"]
                blast_contigs_files_paths.append(ref_genome)
                blast_samples_basenames.append(reference_genome_basename)

            # Blast can run out of memory if there are too many threads and fail without any error message
            blast_state = -1
            blast_n_threads = n_threads
            while blast_n_threads > 0:
                blast_state = blast_call(proteins_file_ori=proteins_file, 
                                        proteins_file_dest=blast_proteins_dir+"/"+proteins_database_name, 
                                        contigs_files_paths=blast_contigs_files_paths, 
                                        blast_database_output=dna_database_blast+"/DNA_database.fna", 
                                        blast_output_folder=blast_proteins_dir,
                                        blast_output_name=blast_output_name,
                                        threads = blast_n_threads)

                if blast_state != 0: 
                    print("ERROR: BLAST failed. Probably because it ran out of memory.", flush=True)
                    print("BLAST exit code: ", blast_state, flush=True)
                    print("Retrying with less threads...", flush=True)
                    blast_n_threads = int(blast_n_threads/2)
                else:
                    break
            if blast_state != 0:
                print("ERROR: BLAST failed.", flush=True)
                sys.exit(1)
                    
            blast_postprocessing(blast_file=blast_proteins_dir+"/"+blast_output_name,
                                database_file=blast_proteins_dir+"/"+proteins_database_name,
                                output_folder=blast_proteins_dir,
                                samples=blast_samples_basenames)

            
    # ABRicate call (Plasmids)
    if cfg.config["plasmids"]["run_plasmid_prediction"]:
        plasmids_database = "plasmidfinder"
        plasmids_output_file = "Plasmids_ABRicate_plasmidfinder.tsv"
        plasmids_matrix_file = "Plasmids_ABRicate_plasmidfinder_matrix.tsv"
        print(Banner(f"\nStep {step_counter}: Plasmids (ABRicate: "+plasmids_database+")\n"), flush=True)
        step_counter += 1
        abricate_samples_basenames = samples_basenames.copy()
        if cfg.config["reference_genome"]["strain"]:
            abricate_samples_basenames.append(cfg.config["reference_genome"]["strain"])
        abricate_call(input_dir=draft_contigs_dir,
                    output_dir=plasmid_dir,
                    output_filename=plasmids_output_file,
                    database=plasmids_database,
                    mincov=cfg.config["plasmids"]["abricate"]["mincov"],
                    minid=cfg.config["plasmids"]["abricate"]["minid"],
                    gene_matrix_file=plasmids_matrix_file,
                    samples=abricate_samples_basenames,
                    threads = n_threads)

        # Delete plasmids_matrix_file if it contains only one line
        with open(plasmid_dir+"/"+plasmids_matrix_file, 'r') as f:
            if len(f.readlines()) == 1:
                print(f"\nINFO: No plasmids were found.\n", flush=True)

        # End line
        with open(plasmid_dir+"/"+plasmids_matrix_file, "a") as matrix_file:
            matrix_file.write("Coverage >= " + 
                                str(cfg.config["plasmids"]["abricate"]["mincov"]) +
                                " % and identity >= " +
                                str(cfg.config["plasmids"]["abricate"]["minid"]) + 
                                " % on each sample for considering a plasmid as present.")


    


    # Antimicrobial resistance genes
    if cfg.config["antimicrobial_resistance_genes"]["run_antimicrobial_resistance_genes_prediction"]:
        if "abricate" in cfg.config["antimicrobial_resistance_genes"]["antimicrobial_resistance_genes_predictor_tool"]:
            global_amr_output_file = "AMR_ABRicate.tsv"
            global_amr_matrix_file = "AMR_genes_ABRicate_matrix.tsv"
            for amr_db in cfg.config["antimicrobial_resistance_genes"]["abricate"]["antimicrobial_resistance_databases"]:
                print(Banner(f"\nStep {step_counter}: Antimicrobial resistance genes (ABRicate: "+amr_db+")\n"), flush=True)
                step_counter += 1
                amr_output_file = "AMR_ABRicate_"+amr_db+".tsv"
                amr_matrix_file = "AMR_genes_ABRicate_"+amr_db+"_matrix.tsv"
                abricate_samples_basenames = samples_basenames.copy()
                if cfg.config["reference_genome"]["strain"]:
                    abricate_samples_basenames.append(cfg.config["reference_genome"]["strain"])
                abricate_call(input_dir=draft_contigs_dir, 
                            output_dir=amr_analysis_dir_abr,
                            output_filename=amr_output_file,
                            database=amr_db,
                            mincov=cfg.config["antimicrobial_resistance_genes"]["abricate"]["mincov"],
                            minid=cfg.config["antimicrobial_resistance_genes"]["abricate"]["minid"],
                            gene_matrix_file=amr_matrix_file,
                            samples=abricate_samples_basenames,
                            threads = n_threads)

                # Concatenate every ABRicate output in a single file
                with open(amr_analysis_dir_abr+"/"+global_amr_output_file, "a") as global_file, open(amr_analysis_dir_abr+"/"+amr_output_file, "r") as current_file:
                    if os.stat(amr_analysis_dir_abr+"/"+global_amr_output_file).st_size == 0:   # If global file is empty, i.e. there is no header
                        global_file.write(current_file.read())
                    else:
                        for line in current_file.readlines()[1:]:
                            global_file.write(line)
                    
                with open(amr_analysis_dir_abr+"/"+global_amr_matrix_file, "a") as global_matrix, open(amr_analysis_dir_abr+"/"+amr_matrix_file, "r") as current_matrix:
                    if os.stat(amr_analysis_dir_abr+"/"+global_amr_matrix_file).st_size == 0:   # If global file is empty, i.e. there is no header
                        global_matrix.write(current_matrix.read())
                    else:
                        for line in current_matrix.readlines()[1:]:
                            global_matrix.write(line)
                
                # Remove innecesary single files
                os.remove(amr_analysis_dir_abr+"/"+amr_output_file)
                os.remove(amr_analysis_dir_abr+"/"+amr_matrix_file)

            # End line
            with open(amr_analysis_dir_abr+"/"+global_amr_matrix_file, "a") as matrix_file:
                matrix_file.write("Coverage >= " + 
                                    str(cfg.config["antimicrobial_resistance_genes"]["abricate"]["mincov"]) +
                                    " % and identity >= " +
                                    str(cfg.config["antimicrobial_resistance_genes"]["abricate"]["minid"]) + 
                                    " % on each sample for considering a virulence gene as present.")
    
        if "amrfinder" in cfg.config["antimicrobial_resistance_genes"]["antimicrobial_resistance_genes_predictor_tool"]:
            amrfinder_db_name = "NDARO"
            print(Banner(f"\nStep {step_counter}: Antimicrobial resistance genes and point mutations (AMRfinder: {amrfinder_db_name})\n"), flush=True)
            step_counter += 1
            if annotator == "prokka":
                gff_dir = prokka_refactor_dir
            elif annotator == "dfast":
                gff_dir = dfast_refactor_dir
            else:
                print("Specified annotator("+annotator+") is not valid.", flush=True)

            amrfinder_out_file = amr_analysis_dir_amrfinder+"/AMR_AMRFinder.tsv"
            amrfinder_resume_file = amr_analysis_dir_amrfinder+"/AMR_genes_AMRFinder_matrix.tsv"
            amrfinder_point_mutations_file = amr_analysis_dir_amrfinder+"/AMR_point_mutations_AMRFinder_matrix.tsv"

            amrfinder_samples = samples_basenames.copy()
            if cfg.config["reference_genome"]["file"]:
                amrfinder_samples.append(reference_genome_basename)
            else:
                reference_genome_basename = None
            amrfinder_call(amrfinder_out_file, amrfinder_resume_file, amrfinder_samples, 
                           annotation_dir, gff_dir, genus, amrfinder_db_name, amr_analysis_dir_amrfinder,
                           ref_genome_basename=reference_genome_basename,
                           threads = n_threads)
            amrfinder_get_point_mutations(amrfinder_out_file, amrfinder_point_mutations_file)

    # Roary call
    if cfg.config["pangenome"]["run_pangenome"]:
        print(Banner(f"\nStep {step_counter}: Roary\n"), flush=True)
        step_counter += 1
        roary_call(input_files=roary_input_files, output_dir=roary_dir, 
                   campype_output_folder=output_folder, threads = n_threads)

        if cfg.config["MLST"]["run_mlst"]:
            # Roary plots call
            os.mkdir(roary_plots_dir)
            print(Banner(f"\nStep {step_counter}: Roary Plots\n"), flush=True)
            step_counter += 1
            roary_plots_call(input_newick=roary_dir+"/accessory_binary_genes.fa.newick",
                            input_gene_presence_absence=roary_dir+"/gene_presence_absence.csv",
                            output_dir=roary_plots_dir)

    # Final report
    print(Banner(f"\nStep {step_counter}: Generate final report\n"), flush=True)
    step_counter += 1
    amrfinder_matrix_file = False
    if cfg.config["antimicrobial_resistance_genes"]["run_antimicrobial_resistance_genes_prediction"]:
            if "amrfinder" in cfg.config["antimicrobial_resistance_genes"]["antimicrobial_resistance_genes_predictor_tool"]:
                amrfinder_matrix_file = amr_analysis_dir_amrfinder+"/AMR_genes_AMRFinder_matrix.tsv"
    
    if cfg.config["summary_report"]["include_reference"] or not cfg.config["reference_genome"]["file"]:
        generate_report(samples_basenames,
                        prinseq_dir,
                        assembly_dir,
                        annotation_dir,
                        draft_contigs_dir,
                        output_folder,
                        summary_pre_qc,
                        summary_post_qc,
                        kraken_unified_report,
                        summary_post_flash,
                        mlst_data_file,
                        blast_proteins_dir+"/Virulence_genes_BLAST_matrix.tsv",
                        snippy_summary=snippy_summary_outfile,
                        custom_VFDB=blast_proteins_dir+"/"+proteins_database_name,
                        amrfinder_matrix_file=amrfinder_matrix_file,
                        reference_genome_basename=reference_genome_basename)
    else:
        generate_report(samples_basenames,                    
                        prinseq_dir, 
                        assembly_dir, 
                        annotation_dir,
                        draft_contigs_dir,
                        output_folder, 
                        summary_pre_qc, 
                        summary_post_qc, 
                        kraken_unified_report,
                        summary_post_flash,
                        mlst_data_file,
                        blast_proteins_dir+"/Virulence_genes_BLAST_matrix.tsv",
                        snippy_summary=snippy_summary_outfile,
                        custom_VFDB=blast_proteins_dir+"/"+proteins_database_name,
                        amrfinder_matrix_file=amrfinder_matrix_file)


    # Remove temporal folders
    if cfg.config["trim_adaptors"]["run_trim_adaptors"]:
        shutil.rmtree(trimmomatic_dir)
    if cfg.config["reference_genome"]["file"]:
        shutil.rmtree(contigs_dir)
    if cfg.config["antimicrobial_resistance_genes"]["run_antimicrobial_resistance_genes_prediction"]:
        if "amrfinder" in cfg.config["antimicrobial_resistance_genes"]["antimicrobial_resistance_genes_predictor_tool"]:
            if annotator == "prokka":
                shutil.rmtree(prokka_refactor_dir)

    # Remove temporal files
    if compressed_mode:
        # Remove decompressed files
        if fasta_mode:
            for key, value in decompressed_samples.items():
                os.remove(value)
        else:
            for key, value in decompressed_samples_fw.items():
                os.remove(value)
            for key, value in decompressed_samples_rv.items():
                os.remove(value)

        # Keep only compressed files in trimmomatic/prinseq folders to reduce the output size
        for root, dirnames, filenames in os.walk(prinseq_dir, topdown=False):
            for dirname in dirnames:
                try:
                    os.remove(os.path.join(root, dirname, dirname+"_R1.fastq"))
                    os.remove(os.path.join(root, dirname, dirname+"_R2.fastq"))
                except OSError:
                    pass
    else:
        # Compressed fastq files in trimmomatic/prinseq folders to reduce the output size
        for root, dirnames, filenames in os.walk(prinseq_dir, topdown=False):
            for dirname in dirnames:
                fastqR1 = os.path.join(root, dirname, dirname+"_R1.fastq")
                fastqR2 = os.path.join(root, dirname, dirname+"_R2.fastq")
                
                for fastqfile in [fastqR1, fastqR2]:
                    if os.path.exists(fastqfile):
                        call(["gzip", fastqfile])
        
    # Gzip flash fastq files to reduce the output size
    for root, dirnames, filenames in os.walk(flash_dir, topdown=False):
        for dirname in dirnames:
            extendedfragsfile = os.path.join(root, dirname, dirname+".extendedFrags.fastq") # TODO should I compress this one too? It was not specified in the mail
            notcombinedfile1 = os.path.join(root, dirname, dirname+".notCombined_1.fastq")
            notcombinedfile2 = os.path.join(root, dirname, dirname+".notCombined_2.fastq")
            
            for fastqfile in [extendedfragsfile, notcombinedfile1, notcombinedfile2]:
                if os.path.exists(fastqfile):
                    call(["gzip", fastqfile])

    # Remove empty folders
    for root, dirnames, filenames in os.walk(output_folder, topdown=False):
        for dirname in dirnames:
            try:
                os.rmdir(os.path.join(root, dirname))
            except OSError:
                pass
