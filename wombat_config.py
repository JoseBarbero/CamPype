config = {
    "output_directory": ".",    # "." is Wombat's directory. (You need to have writing permissions to save Wombat's output in certain directories.)
    "adapters_reference_file": "reference_files/adapters_and_sequences.fa", 
    "reference_genome": {   # OPTIONAL (if you don't want to use it, just leave every parameter empty)
        "file": "reference_files/NCTC11168.fasta", 
        "genus": "Campylobacter",   
        "species": "jejuni",    
        "strain": "NCTC11168",
        "proteins": "reference_files/NCTC11168_NCBI.gb", # TODO optional
    },
    "proteins_reference_file": "reference_files/VF_custom.txt", # OPTIONAL (mandatory if "run_blast" is set to True)
    "run_trimmomatic": True, # Set to True or False,
    "min_contig_len": 500,
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
        "mode": "--isolate",    # Modes: --isolate, --careful or --sc. Default: --isolate
        "cov_cutoff": "auto",    # Read coverage cutoff value. Must be a positive float value, or 'auto', or 'off'. Default value is 'auto', when SPAdes automatically computes coverage threshold using conservative strategy
        "k": False     # k-mer sizes to be used. Set this to a number or comma-separated list of numbers (all values must be odd, less than 128 and listed in ascending order). Example: 55,77,121. If set to False, values are automatically selected using maximum read length. If mode --sc is set, the default values are 21, 33 and 55. Default: False 
    },
    "quast": {
        "icarus": "--no-icarus",    # Do not build Icarusviewers
        "mode": "--silent"  # Do not print detailed information about each step in standard output. This option does not affect quast.log file
    },
    "annotator": "prokka",  # Set this to "prokka" or "dfast"
    "prokka": {
        "kingdom": "Bacteria",  # Annotation mode: Archaea|Bacteria|Mitochondria|Viruses (default 'Bacteria')
        "gcode": 11,  # Genetic code / Translation table (set if --kingdom is set) (default '11')
        "metagenome": True, # Set to True or False TODO 
        "rawproduct": True, # Set to True or False TODO 
        "reference_annotation": True    # Use reference genome annotation GenBank file to first annotate from. Default: True
    },
    "dfast": {
        "min_length": 1,    # Minimum sequence length. Default: 1
        "use_original_name": "true",    # Use original sequence names in a query FASTA file- Default: True
        "sort": "false",
        "step": 1
    },
    "abricate": {
        "run_amr": True,  # Run antimicrobial resistance gene search. If you want to omit this step, set it to False. Default: False
        "virus_database": "vfdb",
        "antimicrobial_resistance_database": "card", # Select database from ABRicate to identify antimicrobial resistance genes: argannot, card, ecoh, ecoli_vf, megares, ncbi or resfinder. Remember, ABRicate uses blastn. Default: card
        "mincov": 90, # Minimum DNA % coverage for considering an antimicrobial resistance gene as present. Default: 90
        "minid": 90 # Minimum DNA % identity for considering an antimicrobial resistance gene as present. Default: 90      
    },
    "amrfinder": {
        "run": True,  # If you want to omit this step, set it to False. Default: True
        "update_db": False,  # Updates the db before running amrfinder. Set to True or False. Default: False
        "minid": 0.9,       # Minimum proportion identical translated AA residues for considering an antimicrobial resistance gene (0-1). Default: 0.9
        "mincov": 0.9       # Minimum coverage of reference protein sequence for for considering an antimicrobial resistance gene (0-1). Default: 0.9
    },
    "run_blast": True,     # Set this to True or False
    "blast": {
        "dbtype": "nucl",   # TODO
        "evalue": 10e-4,    # TODO
        "outfmt": "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq", # Output format
        "soft_masking": True    # Apply filtering locations as soft masks (i.e., only for finding initial matches)
    },
    "presence_absence_matrix": {
        "protein_cover": 90,    # Minimum % protein cover for considering a virulence gene as present. Default: 90
        "protein_identity": 90  # Minimum % protein identity for considering a virulence gene as present. Default: 90
    },
    "roary":{
        "split_paralogs": True, # Set this to True (don't split paralogs) or False (split paralogs). Default: True
        "min_identity": 95  # Minimum percentage identity for blastp (1-100). Default: 95
    }

}
