config = {
    "output_directory": ".",    # "." is Wombat's directory. (You need to have writing permissions to save Wombat's output in certain directories.)
    "adapters_reference_file": "reference_files/adapters.fa", # TODO make this generic for other users
    "annotation_reference": {
        "file": "reference_files/NCTC11168.fasta", # TODO make this generic for other users
        "genus": "Campylobacter",   # TODO make this generic for other users
        "species": "jejuni",    # TODO make this generic for other users
        "strain": "NCTC11168" # TODO make this generic for other users
    },
    "proteins_reference_file": "reference_files/VF_custom.txt", # TODO make this generic for other users
    "run_trimmomatic": False, # Set to true or false,
    "min_contig_len": 500,
    "prinseq": {
        "min_len": 40, # Minimum read length. (default: {40})
        "min_qual_mean": 25, # Minimum read quality. (default: {25})
        "trim_qual_right": 25,  # Trim sequence by quality score from the 3'-end with this threshold score. (default: {25})
        "trim_qual_window": 15,     # The sliding window size used to calculate quality score by type. To stop at the first base that fails the rule defined, use a window size of 1
        "trim_qual_type": "mean",   # Type of quality score calculation to use. (default: {"mean"})
        "out_format": 3,    # Output format 1 (FASTA only), 2 (FASTA and QUAL), 3 (FASTQ), 4 (FASTQ and FASTA), 5 (FASTQ, FASTA and QUAL) (default: {3})
        "out_bad": "null"   # Data not passing anyfilter will be ignored.
        },  
    "spades": {
        "mode": "--careful",    # Minimize number ofmismatches in the final contigs.
        "cov_cutoff": "auto"    # Read coverage cutoff value. Must be a positive float value, or 'auto', or 'off'. Default value is 'auto', when SPAdes automatically computes coverage threshold using conservative strategy. 
    },
    "quast": {
        "icarus": "--no-icarus",    # Do not build Icarusviewers.
        "mode": "--silent"  # Do not print detailed information about each step in standard output. This option does not affect quast.log file.
    },
    "annotator": "dfast",  # Set this to "prokka" or "dfast"
    "prokka": {
        "kingdom": "Bacteria",  # Annotation mode: Archaea|Bacteria|Mitochondria|Viruses (default 'Bacteria')
        "gcode": 11 # Genetic code / Translation table (set if --kingdom is set) (default '11')
    },
    "dfast": {
        "min_length": 0,    # Minimum sequence length (default '0').
        "use_original_name": "true",    # Use original sequence names in a query FASTA file (default 'true')
        "sort": "false",
        "step": 1
    },
    "mlst": {},
    "abricate": {
        "virus_database": "vfdb",
        "bacteria_database": "resfinder"
    },
    "run_blast": True,     # Set this to True or False
    "blast": {
        "dbtype": "nucl",   # TODO
        "evalue": 10e-4,    # TODO
        "outfmt": "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" # Output format
    }
}