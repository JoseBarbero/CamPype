# Wombat

Wombat is a pipeline for the analysis of Illumina paired-end sequencing data and/or whole bacterial genomes. The development of the workflow is mainly intended for the analysis of Campylobacter jejuni/coli genomes, although any other bacterial genus can be analyzed as well. Wombat is specially designed for users without knowledge of bioinformatics or programming, so the ease of installation and execution are the fundamentals of its development. Moreover, Wombat is a user-customizable workflow that allows you to select the analysis and the tools you are interested in.

![scheme](https://user-images.githubusercontent.com/58036036/180643838-d771f326-3ef9-465e-b591-b5e9df792aec.png)

## Installation (Linux)

1. Clone this repository
    ```bash
    $ git clone https://github.com/JoseBarbero/Wombat.git
    ```
1. Go to Wombat's folder
    ```bash
    $ cd Wombat
    ```
1. Create the environment with conda
    ```bash
    $ conda config --append channels conda-forge
    $ conda config --append channels bioconda
    $ conda env create -f wombatenv.yml 
    $ conda env create -f wombatenv_aux.yml 
    ```
1. *(On Ubuntu systems you may find missing font problems running Mauve. We recommend you to install the required fonts just in case.)
    ```bash
    $ sudo apt-get install ttf-dejavu
    ```

IMPORTANT. After installing or updating Wombat, we recommend you to update the databases of AMRFinder, Prokka and ABRicate by running:
 ```
 conda activate wombat
 amrfinder -u
 prokka --setupdb
 abricate --setupdb
 ```
 
## How to update Wombat (Linux)

1. Make sure you are not in Wombat conda environment
    ```bash
    $ conda deactivate
    ```
1. Save a copy of your configuration files. Updating Wombat will overwrite you configuration files because this files properties may change with newer versions of Wombat.

1. Run ./updatewombat
    ```bash
    $ ./updatewombat
    ```
We recommend you to answer YES to the first question to update configuration files and NO to the second question, otherwise all the packages and tools included in Wombat will me removed and installed again, which is a time-consuming task.


## How to uninstall Wombat (Linux)

1. Make sure you are not in Wombat conda environment
    ```bash
    $ conda deactivate
    ```
1. Every file in Wombat folder will be deleted. Make sure you don't have anything important in this directory (results, data, etc).

1. Run ./uninstallwombat
    ```bash
    $ ./uninstallwombat
    ```


## Set input files and configuration

1.  Set input files in Wombat/input_files.csv (tab as separator):

    | Samples        | Forward           | Reverse  | Genus  | Species  |
    | ------------- |:-------------:|:-----:|:-----:|:-----:|
    | [sample_basename]  | /path/to/your/forward/fastq1_file.fastq | /path/to/your/reverse/fastq1_file.fastq | YourStrainGenus | YourStrainSpecies
    | [sample_basename]  | /path/to/your/forward/fastq2_file.fastq | /path/to/your/reverse/fastq2_file.fastq | YourStrainGenus | YourStrainSpecies
    | [sample_basename]  | /path/to/your/forward/fastq3_file.fastq | /path/to/your/reverse/fastq3_file.fastq | YourStrainGenus | YourStrainSpecies

    Make sure headers haven't changed. If your input data are previously assembled genomes, both the Forward and Reverse columns will refer to the path of the assembled genomes, that is, the content of both columns will be the same, except for the headers name.
    
1. Set your own running parameters and tools to run in "workflow_config.py" file. 

1. Default settings are configured for Campylobacter jejuni/coli. I you want to use a different bacteria, we recommend you to deactivate the option ```include_cc``` and use abricate for virulence genes searching or use your own database with BLAST instead.


## Running the workflow

1. Activate the environment
    ```bash
    $ conda activate wombat
    ```
1. Go to Wombat folder
    ```bash
    $ cd your/path/to/Wombat
    ```
1. Run the workflow
    ```bash
    $ bash -i wombat
    ```
1. \(*) You can deactivate the environment when you are finished
    ```bash
    $ conda deactivate
    ```

## Output
The results of Wombat are stored in very detailed directories for each analysis, with separate folders for each tool and isolate. Los files will be generated for analysis tracking due to execution error. An interactive HTML summary report will be generated at the end of the analysis to simplify the task of data visualization and interpretation. This HTML file can be opened on any Web browser. An example of report can be found here.

You can generate the report after Wombat analysis by executing the following command in the Linux terminal:
```
Rscript -e "rmarkdown::render('Wombat_Report_long.Rmd', params = list(directory = '~/path/to/data'))"
Rscript -e "rmarkdown::render('Wombat_Report_short.Rmd', params = list(directory = '~/path/to/data'))"
```
where you will have to change ```'~/path/to/data'``` with the corresponding path of the Wombat output folder.


## FAQ
1. Prokka stops running with this error:
```
Could not run command: cat \/home\/Wombat_OUTPUT_20220511_131550\/Prokka_annotation\/NCTC11168\/NCTC11168\.IS\.tmp\.35844\.faa | parallel --gnu --plain -j 8 --block 313 --recstart '>' --pipe blastp -query - -db /home/instalador/anaconda3/envs/wombat/db/kingdom/Bacteria/IS -evalue 1e-30 -qcov_hsp_perc 90 -num_threads 1 -num_descriptions 1 -num_alignments 1 -seg no > \/home\/Wombat_OUTPUT_20220511_131550\/Prokka_annotation\/NCTC11168\/NCTC11168\.IS\.tmp\.35844\.blast 2> /dev/null
```

Run ```prokka --setupdb``` first and execute Wombat again.

2.  ABRicate can't find any gen and this message appears: ```BLAST Database error: Error pre-fetching sequence data```

Run ```abricate --setupdb``` first and execute Wombat again.


## Citation
Please cite Wombat whenever you use it as:

Irene Ortega-Sanz, Jose A. Barbero and Antonio Canepa. Wombat. Availabe at [https://github.com/JoseBarbero/Wombat](https://github.com/JoseBarbero/Wombat)


## Contact:
For questions, bugs and suggestions, please contact to iortega@ubu.es or jabarbero@ubu.es.

