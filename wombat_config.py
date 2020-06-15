config = {
    "output_directory": ".",    # "." is Wombat's directory (you need to have writing permissions to save Wombat's output in certain directories.)
    "adapters_reference_file": "reference_files/adapters_and_sequences.fa",   # You can edit this file with any sequence you want to filter
    "reference_genome": {   # OPTIONAL (if you don't want to use it, just leave every parameter empty)
        "file": "reference_files/NCTC11168.fasta",    # Required if you want to reorder genome contigs against a reference genome and search SNPs
        "genus": "Campylobacter",   
        "species": "jejuni",    
        "strain": "NCTC11168",
        "proteins": "reference_files/NCTC11168_NCBI.gb", # Required if you want to reduce annotation mismatches using PROKKA compared to a reference genome and get detailed information of SNPS 
    },
    "run_trimmomatic": True, # Set to True or False if you want to run Trimmomatic or not to filter sequenced reads
    "prinseq": {
        "min_len": 50, # Minimum read length. Default: 50
        "min_qual_mean": 30, # Minimum read quality. Default: 30
        "trim_qual_right": 25,  # Trim sequence by quality score from the 3'-end with this threshold score. Default: 25
        "trim_qual_window": 20,     # The sliding window size used to calculate quality score by type. To stop at the first base that fails the rule defined, use a window size of 1. Default: 20
        "trim_qual_type": "mean",   # Type of quality score calculation to use. Default: 20
        "out_format": 3,    # Output format 1 (FASTA only), 2 (FASTA and QUAL), 3 (FASTQ), 4 (FASTQ and FASTA), 5 (FASTQ, FASTA and QUAL) (default: {3})
        "out_bad": "null"   # Data not passing anyfilter will be ignored
    },  
    "spades": {
        "mode": "--isolate",    # Assembly modes: --isolate (highly recommended for high-coverage isolate and multi-cell Illumina data), --careful (reduce the number of mismatches and short indels, recommended only for assembly of small genomes) or --sc (required for MDA (single-cell) data). Default: --isolate
        "cov_cutoff": "auto",    # Read coverage cutoff value. Must be a positive float value, or 'auto', or 'off'. Default value is 'auto', when SPAdes automatically computes coverage threshold using conservative strategy
        "k": False     # k-mer sizes to be used. Set this to a number or comma-separated list of numbers (all values must be odd, less than 128 and listed in ascending order). Example: 55,77,121. If set to False, values are automatically selected using maximum read length. If mode --sc is set, the default values are 21, 33 and 55. Default: False 
    },
    "min_contig_len": 500,
    "quast": {
        "icarus": "--no-icarus",    # Do not build Icarusviewers
        "mode": "--silent"  # Do not print detailed information about each step in standard output. This option does not affect quast.log file
    },
    "annotator": "prokka",  # Set this to "prokka" or "dfast"
    "prokka": {
        "kingdom": "Bacteria",  # Annotation mode: Archaea|Bacteria|Mitochondria|Viruses (default 'Bacteria')
        "gcode": 11,  # Genetic code / Translation table (set if --kingdom is set) (default '11')
        "metagenome": True, # Improve gene predictions for highly fragmented genomes if set to True (not for reference genome). Otherwise, set to False. Default: True
        "rawproduct": True, # Do not clean up product annotation if set to True. Otherwise, set to False. Default: True
        "reference_annotation": True    # Use reference genome annotation GenBank file to first annotate from if set to True. Otherwise, set to False. Default: True
    },
    "dfast": {
        "min_length": 1,    # Minimum sequence length (bp) to annotate. Default: 1
        "use_original_name": "true",    # Use original sequence names in a query FASTA file (true or false). Default: True
        "sort": "false",    # Sort sequences by length (true or false). Default: False
        "step": 1
    },
    "abricate": {
        "run_amr": True,  # Run antimicrobial resistance gene search using ABRicate. If you want to omit this step, set it to False. Default: False
        "virus_database": "vfdb",
        "antimicrobial_resistance_database": "card", # Select database from ABRicate to identify antimicrobial resistance genes: argannot, card, ecoh, ecoli_vf, megares, ncbi or resfinder. Remember, ABRicate uses blastn. Default: card
        "mincov": 90, # Minimum DNA % coverage for considering an antimicrobial resistance gene as present. Default: 90
        "minid": 90 # Minimum DNA % identity for considering an antimicrobial resistance gene as present. Default: 90      
    },
    "amrfinder": {
        "run": True,  # Run antimicrobial resistance gene search using AMRFinder (database: NDARO). If you want to omit this step, set it to False. Default: True
        "update_db": False,  # Updates the database before running AMRFinder (Internet connection is required) if set to True. Otherwise, set to False. Default: False
        "minid": 0.9,       # Minimum proportion identical translated AA residues for considering an antimicrobial resistance gene (0-1). Default: 0.9
        "mincov": 0.9       # Minimum coverage of reference protein sequence for for considering an antimicrobial resistance gene (0-1). Default: 0.9
    },
    "run_blast": True,     # Run tBLASTn to scan specific virulence genes from the custom_VFDB.txt. Remember that you can edit that database with your own sequences. If you want to omit this step, set it to False. Default: True
    "proteins_reference_file": "reference_files/VF_custom.txt", # OPTIONAL (mandatory if "run_blast" is set to True). You can edit this file with any sequence you want to scan
    "blast": {
        "dbtype": "nucl",   # TODO
        "evalue": 10e-4,    # Maximum evalue for a hit. Default: 10e-4
        "outfmt": "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq", # Output format
        "soft_masking": True    # Apply filtering locations as soft masks (i.e., only for finding initial matches) if set to True. Otherwise, set to False. Default: True
    },
    "presence_absence_matrix": {    # Parameteres for blast matrix construction
        "protein_cover": 90,    # Minimum % protein cover for considering a virulence gene as present. Default: 90
        "protein_identity": 90  # Minimum % protein identity for considering a virulence gene as present. Default: 90
    },
    "roary":{
        "split_paralogs": True, # Set this to True (do not split paralogs) or False (split paralogs). Default: True
        "min_identity": 95  # Minimum percentage identity for blastp (1-100). Default: 95
    }

}
