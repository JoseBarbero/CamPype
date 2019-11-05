# Wombat workflow

## Installation (Linux)

1. Clone this repository
    ```bash
    $ git clone https://github.com/JoseBarbero/Wombat.git
    ```
1. Go to Wombat's folder
    ```bash
    $ cd Wombat
    ```
1. Create a conda environment
    ```bash
    $ conda env create -f wombatenv.yml
    ```
    
## Set input files 

1.  Set input files in Wombat/input_files.csv (tab as separator):

    | Samples        | Read1           | Read2  |
    | ------------- |:-------------:| -----:|
    | [sample_basename]  | /path/to/your/forward/fastq1_file.fastq | /path/to/your/reverse/fastq1_file.fastq |
    | [sample_basename]  | /path/to/your/forward/fastq2_file.fastq | /path/to/your/reverse/fastq2_file.fastq |
    | [sample_basename]  | /path/to/your/forward/fastq3_file.fastq | /path/to/your/reverse/fastq3_file.fastq |
   
1.  Set auxiliary files in Wombat/reference_files.csv (tab as separator):

    | File        | Route           |
    | ------------- |:-------------:|
    | adapters  | path/to/your/sequencing/adapters/file/adapters.fa |
    | reference_annotation_file  | path/to/your/reference/genome/file/NCTC11168.fasta |

## Running the workflow

1. Activate the environment
    ```bash
    $ conda activate wombat
    ```
1. Go to Wombat's folder
    ```bash
    $ cd your/path/to/Wombat
    ```
1. Run the workflow
    ```bash
    $ python wombat.py
    ```
1. \(*) You can deactivate the environment when you are finished
    ```bash
    $ conda deactivate
    ```
