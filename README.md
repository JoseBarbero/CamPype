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
1. *(On Ubuntu systems you may find missing font problems running Mauve. We recomend you to install the required fonts just in case.)
    ```bash
    $ sudo apt-get install ttf-dejavu
    ```

## How to update Wombat (Linux)

1. Make sure you are not inside Wombat conda environment
    ```bash
    $ conda deactivate
    ```
1. Save a copy of your configuration files. Updating Wombat will overwrite you configuration files because this files properties may change with newer versions.

1. Run ./updatewombat
    ```bash
    $ ./updatewombat
    ```

## How to uninstall Wombat (Linux)

1. Make sure you are not inside Wombat conda environment
    ```bash
    $ conda deactivate
    ```
1. Every file into Wombat folder will be deleted. Make sure you don't have anything important in this directory.

1. Run ./uninstallwombat
    ```bash
    $ ./uninstallwombat
    ```





## Set input files and configuration

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
    
    \* Remove the last line if you don't have a reference file to annotate.

1. Set your own running parameters in "workflow_config.py" file.

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
    $ ./wombat
    ```
1. \(*) You can deactivate the environment when you are finished
    ```bash
    $ conda deactivate
    ```
