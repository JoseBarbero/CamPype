config = {
    "output_directory": ".",    # "." is CamPype's directory (you need to have writing permissions to save CamPype's output in certain directories.). Default: "."
    "custom_output_name": "", # Custom name for the output directory. Default: "" (if this parameter is empty: "CamPype_OUTPUT_YYYYMMDD_HHMMSS")
    "n_threads": 28,   # Number of threads to use. Default: None (takes threads defined by default by each tool). Default: 28
    "reference_genome": {   # OPTIONAL (if you don't want to use it, just leave every parameter empty)
        "file": "reference_files/NCTC11168.fasta",    # Required if you want to reorder genome contigs against a reference genome and search for SNPs
        "genus": "Campylobacter",   
        "species": "jejuni",    
        "strain": "NCTC11168",
        "proteins": "reference_files/NCTC11168_NCBI.gb", # Required if you want to reduce annotation mismatches using PROKKA compared to a reference genome and get detailed information of SNPS. Genome sequence is required at the end of the file for successful execution 
    },
    "assembled_genomes": False, # Skip read quality control and assembly if set to True (you need fasta files as input). If not, set it to False. Default: False
    "species_identification": { 
        "run_species_identification": False, # Skip bacteria classification if you already know the species you are working with. In case, you set it to True and you already indicated the Genus and Species of your genomes in the input_files.csv, the output of this step will not be used in further analysis. Default: False
        "species_identification_database": "./db/minikraken-DB/minikraken_8GB_20200312" # Location of standardized database for bacteria classification. Default: "./db/minikraken-DB/minikraken_8GB_20200312"
    },
    "trim_adaptors":{
      "run_trim_adaptors": True, # Remove adapter sequences from reads if set to True. Otherwise, set to False. Default: True
      "adapters_reference_file": "reference_files/adapters_and_sequences.fa"   # Location of the sequences to filter. You can edit this file with any sequence you want to filter. Default: "reference_files/adapters_and_sequences.fa"
    },
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
    "min_contig_len": 200,
    "MLST": {
        "run_mlst": True, # Set this to True or False if you want to run the MLST step. Default: True
        "include_cc": True # Will assign Sequence Types (ST) to Clonal Complexes (CC) using the PubMLST Campylobacter jejuni/coli database (internet connection is required). You must set this option to False if you are analysing other genus. Default: True
    },
    "annotation": {
        "run_annotation": True, # Set this to True or False if you want to run the annotation step. Default: True
        "annotator": "prokka",  # Set this to "prokka" or "dfast". Default: prokka
        "prokka": {
            "rawproduct": True, # Do not clean up product annotation if set to True. Otherwise, set to False. Default: True
            "reference_annotation": True    # Use reference genome annotation GenBank file to first annotate from if set to True. Otherwise, set to False. Default: True
        },
    },
    "antimicrobial_resistance_genes":{
        "run_antimicrobial_resistance_genes_prediction": True, #Set this to True or False if you want to search for antimicrobial resistance genes or not. Default: True
        "antimicrobial_resistance_genes_predictor_tool": ["abricate", "amrfinder"], # Select as many tools as desired to identify virulence genes: abricate and/or amrfinder. Keep in mind that abricate uses BLASTn and AMRFinder uses blastp. Be aware to run genome annotation if using AMRFinder! For Campylobacter jejuni and Campylobacter coli strains, resistance caused by point mutations will be screened as well. Default: ["abricate", "amrfinder"] 
        "abricate": {
            "antimicrobial_resistance_databases": ["argannot", "card", "megares", "ncbi", "resfinder"], # Select as many databases as desired from ABRicate to identify antimicrobial resistance genes: argannot, card, ecoh, ecoli_vf, megares, ncbi or/and resfinder. Remember, ABRicate uses blastn. Default: ["argannot", "card", "megares", "ncbi", "resfinder"]
            "mincov": 80, # Minimum DNA % coverage for considering a virulence gene as present. Default: 80
            "minid": 60 # Minimum DNA % identity for considering a virulence gene as present. Default: 60
        },
        "amrfinder": {
            "mincov": 80, # Minimum coverage of reference protein sequence for for considering an antimicrobial resistance gene (0-100). Default: 80
            "minid": 60 # Minimum proportion identical translated AA residues for considering an antimicrobial resistance gene (0-100). Default: 60
        }
    },
    "virulence_genes":{
        "run_virulence_genes_prediction": True, #Set this to True or False if you want to search for virulence genes or not. Default: True
        "virulence_genes_predictor_tool": ["abricate", "blast"], # Select as many tools as desired to identify virulence genes: abricate and/or blast. Keep in mind that abricate uses BLASTn and our BLAST against an inhouse database uses tBLASTn. Default: ["abricate", "blast"].
        "abricate": {
            "virulence_factors_databases": ["vfdb"], # Select as many databases as desired from ABRicate to identify virulence factors: vfdb or/and ecoli_vf. Remember, ABRicate uses blastn. Default: ["vfdb"]
            "mincov": 80, # Minimum DNA % coverage for considering a virulence gene as present. Default: 80
            "minid": 60 # Minimum DNA % identity for considering a virulence gene as present. Default: 60
        },
        "blast": {
            "proteins_reference_file": "reference_files/Campylobacter_custom_VFDB.txt", # A FASTA database must be indicated if "predictor_tool" above is set to blast. You can modify this inhouse Campylobacter spp. database with any selected protein sequence you want to scan or create a new database for other species. Default: "reference_files/Campylobacter_custom_VFDB.txt" 
            "soft_masking": True, # Apply filtering locations as soft masks (i.e., only for finding initial matches) if set to True. Otherwise, set to False. Default: True
            "presence_absence_matrix": { # Parameters for blast matrix construction
                "mincov": 80, # Minimum % protein cover for considering a virulence gene as present. Default: 80
                "minid": 60 # Minimum % protein identity for considering a virulence gene as present. Default: 60
            }
        }
    },
    "plasmids":{
        "run_plasmid_prediction": True, #Set this to True or False if you want to search for plasmids or not. Default: True
        "abricate": {
            "mincov": 80, # Minimum DNA % coverage for considering a plasmid as present. Default: 80
            "minid": 60 # Minimum DNA % identity for considering a plasmid as present. Default: 60
        }
    },
    "run_variant_calling": True, # Set this to True or False if you want to run the variant calling step. Default: True
    "pangenome":{
        "run_pangenome": True, # Set this to True or False if you want to run the pangenome construction step. Remember that "run_annotation" must be to True to allow the pangenome construction. Default: True
        "no_split_paralogs": True, # Set this to True (do not split paralogs) or False (split paralogs). Default: True
        "minid": 95  # Minimum percentage identity for blastp (1-100). Default: 95
    },
    "summary_report":{
        "include_reference": True, # Set this to True to include the reference genome in the summary report. Default: True
    }
}
