config = {
    "output_directory": ".",    # "." is Wombat's directory (you need to have writing permissions to save Wombat's output in certain directories.)
    "adapters_reference_file": "reference_files/adapters_and_sequences.fa",   # You can edit this file with any sequence you want to filter
    "reference_genome": {   # OPTIONAL (if you don't want to use it, just leave every parameter empty)
        "file": "reference_files/NCTC11168.fasta",    # Required if you want to reorder genome contigs against a reference genome and search SNPs
        "genus": "Campylobacter",   
        "species": "jejuni",    
        "strain": "NCTC11168",
        "proteins": "reference_files/NCTC11168_NCBI.gb", # Required if you want to reduce annotation mismatches using PROKKA compared to a reference genome and get detailed information of SNPS. Genome sequence is required at the end of the file for successful execution 
    },
    "assembled_genomes": True, # Skip read quality control and assembly if set to True (you need fasta files as input). If not, set it to False. Default: False
    "trim_adaptors": False, #Remove adapter sequences from reads if set to True. Otherwise, set to False. Default: True
    "read_qc_filtering": {
        "min_len": 50, # Minimum read length. Default: 50
        "min_qual_mean": 30, # Minimum read quality. Default: 30
        "trim_qual_right": 25,  # Trim sequence by quality score from the 3'-end with this threshold score. Default: 25
        "trim_qual_window": 20     # The sliding window size used to calculate quality score by type. To stop at the first base that fails the rule defined, use a window size of 1. Default: 20
    },
    "merge_reads": True, # Merge pairs of reads before assembly to improve the quality of assembly when paired reads overlap if set to True. Otherwise, set to False. Default: True  
    "assembly": {
        "mode": "--isolate",    # Assembly modes: --isolate (highly recommended for high-coverage isolate and multi-cell Illumina data), --careful (reduce the number of mismatches and short indels, recommended only for assembly of small genomes) or --sc (required for MDA (single-cell) data). Default: --isolate
        "k": False     # k-mer sizes to be used. Set this to a number or comma-separated list of numbers (all values must be odd, less than 128 and listed in ascending order). Example: 55,77,121. If set to False, values are automatically selected using maximum read length. If mode --sc is set, the default values are 21, 33 and 55. Default: False 
    },
    "min_contig_len": 500,
    "MLST": {
        "run_mlst": False, # Set this to True or False if you want to run the MLST step
        "include_cc": False # Will match sequences types (ST) will clonal complexes (CC) using the PubMLST Campylobacter database. You should set this option to False if you are analysing other genus.
    },
    "abricate": {
        "run_amr": True,  # Run antimicrobial resistance gene search using ABRicate and blastn. If you want to omit this step, set it to False. Default: False
        "antimicrobial_resistance_databases": ["card", "ncbi"], # Select as many databases as desired from ABRicate to identify antimicrobial resistance genes: argannot, card, ecoh, ecoli_vf, megares, ncbi and resfinder. Remember, ABRicate uses blastn. Default: card
        "virulence_factors_databases": ["vfdb", "ecoli_vf"], # Select as many databases as desired from ABRicate to identify virulence factors: vfdb and ecoli_vf. Remember, ABRicate uses blastn. Default: vfdb
        "find_plasmids": True,  # Runs ABRicate using plasmidfinder db
        "mincov": 90, # Minimum DNA % coverage for considering an antimicrobial resistance gene as present. Default: 90
        "minid": 90 # Minimum DNA % identity for considering an antimicrobial resistance gene as present. Default: 90      
    },
    "amrfinder": {
        "run": True,  # Run antimicrobial resistance gene search using AMRFinder and blastp (database: NDARO). If you want to omit this step, set it to False. Default: True
        "minid": 0.9,       # Minimum proportion identical translated AA residues for considering an antimicrobial resistance gene (0-1). Default: 0.9
        "mincov": 0.9       # Minimum coverage of reference protein sequence for for considering an antimicrobial resistance gene (0-1). Default: 0.9
    },
    "blast": {
        "run_blast": True,     # Run tBLASTn to scan specific virulence genes from the custom_VFDB.txt. Remember that you can edit that database with your own sequences. If you want to omit this step, set it to False. Default: True
        "proteins_reference_file": "reference_files/custom_VFDB.txt", # OPTIONAL (mandatory if "run_blast" is set to True). You can edit this file with any sequence you want to scan
        "soft_masking": True,    # Apply filtering locations as soft masks (i.e., only for finding initial matches) if set to True. Otherwise, set to False. Default: True
        "presence_absence_matrix": {    # Parameters for blast matrix construction
            "protein_cover": 90,    # Minimum % protein cover for considering a virulence gene as present. Default: 90
            "protein_identity": 90  # Minimum % protein identity for considering a virulence gene as present. Default: 90
        },
    },
    "annotation": {
        "run_annotation": True, # Set this to True or False if you want to run the annotation step
        "annotator": "prokka",  # Set this to "prokka" or "dfast"
        "prokka": {
            "rawproduct": True, # Do not clean up product annotation if set to True. Otherwise, set to False. Default: True
            "reference_annotation": True    # Use reference genome annotation GenBank file to first annotate from if set to True. Otherwise, set to False. Default: True
    },
    "pangenome":{
        "run_pangenome": True, # Set this to True or False if you want to run the pangenome construction step. Remember that "run_annotation" must be to True to allow the pangenome construction
        "split_paralogs": True, # Set this to True (do not split paralogs) or False (split paralogs). Default: True
        "min_identity": 95  # Minimum percentage identity for blastp (1-100). Default: 95
    },
    "report":{
        "include_reference": True, # Set this to True to include the reference genome in the final report.
    }

}
