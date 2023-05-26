# CamPype

CamPype is a pipeline for the analysis of Illumina paired-end sequencing data and/or whole bacterial genomes. The development of the workflow is mainly intended for the analysis of <em>Campylobacter jejuni/coli</em> genomes, although any other bacterial genus can be analyzed as well. CamPype is specially designed for users without knowledge of bioinformatics or programming, so the ease of installation and execution are the fundamentals of its development. Moreover, CamPype is a user-customizable workflow that allows you to select the analysis and the tools you are interested in. 
Here you wil find the schema of CamPype. Software or databases are indicated in boxes, while discontinuous boxes indicate tools that users can deactivate.

![scheme](https://user-images.githubusercontent.com/58036036/180643838-d771f326-3ef9-465e-b591-b5e9df792aec.png)

## Installation (Linux)

1. Clone this repository:
    ```bash
    $ git clone https://github.com/JoseBarbero/CamPype.git
    ```
1. Go to CamPype's directory:
    ```bash
    $ cd CamPype
    ```
1. Create the environment with conda:
    ```bash
    $ conda config --append channels conda-forge
    $ conda config --append channels bioconda
    $ conda env create -f campype_env.yml 
    $ conda env create -f campype_env_aux.yml 
    ```
    The creation of the conda environments will take some minutes, be patient.
    
1. *(On Ubuntu systems you may find missing font problems running Mauve. We recommend you to install the required fonts just in case):
    ```bash
    $ sudo apt-get install ttf-dejavu
    ```
1. Additionally, CamPype allows you to check for read contamination and determine bacteria taxonomy using Kraken2. The installation of Kraken2 is optional to avoid possible storage limitations as it requires the use of a heavy database that requires high free disk space, but it is not needed for Campype if you are not interested in this analysis. If you want to install this module, at least 8 GB will be occupied to store the MiniKraken_8GB_202003 database. This database is enough for bacteria identification and CamPype performance, but if you have enough disk space, you can download and install "Standard Kraken2 Databases" following the instructions [here](https://lomanlab.github.io/mockcommunity/mc_databases.html) for better sensitivity. To install Kraken2 in CamPype run:
    ```bash
    $ install_kraken.sh
    ```

#### $\textcolor{red}{\textsf{IMPORTANT!!}}$ After installing or updating CamPype, we recommend you to update the databases of AMRFinder, Prokka and ABRicate by running:
```
$ conda activate campype
$ amrfinder -u
$ prokka --setupdb
$ abricate --setupdb
```
 
## Set input files and configuration

CamPype can run on two modes depending on the input files. The FASTQ mode analyses (un)compressed raw reads in fastq format, while the FASTA mode analyses assembled genomes in fasta format.

1.  Before running CamPype, you must indicate the location of the input files in the CamPype/input_files.csv file (TAB as separator). For fastq files, you need to indicate the path of each pair of reads in the Forward and Reverse columns, while for fasta files both the Forward and Reverse columns will refer to the path of the assembled genome, that is, the content of both columns will be the same. We recommend you to indicate the Genus and Species of your genomes in case you known this information beforehand (this taxonomy will be considered even in you run the bacteria identification module). Please, indicate the species in case you indicate the genus. If not, we strongly encourage you to activate the bacteria identification analysis (as below explained).

    | Samples        | Forward           | Reverse  | Genus  | Species  |
    | ------------- |:-------------:|:-----:|:-----:|:-----:|
    | sample_ID  | /path/to/your/forward/fastq1_file.fastq | /path/to/your/reverse/fastq1_file.fastq | YourStrainGenus | YourStrainSpecies
    | sample_ID  | /path/to/your/forward/fastq2_file.fastq | /path/to/your/reverse/fastq2_file.fastq | YourStrainGenus | YourStrainSpecies
    | sample_ID  | /path/to/your/forward/fastq3_file.fastq | /path/to/your/reverse/fastq3_file.fastq | YourStrainGenus | YourStrainSpecies

    This structure must be respected anyway in that file. Make sure headers haven't changed! $\textcolor{red}{\textsf{Be careful with typos!!!}}$
    
1. Set the modules you want to run in the CamPype/workflow_config.py file. There you will set your own running parameters for each tool and the select your tools of interest when possible.

1. Default settings are configured for <em>Campylobacter jejuni/coli</em>. If you want to use a different bacteria, we strongly recommend you to adapt the configuration of CamPype as previoys explained. In particular, you must modify the ```reference_genome```, deactivate the option ```include_cc```, use abricate for virulence genes searching or/and use your own virulence genes database with BLAST instead (indicate this accordingly in the CamPype/workflow_config.py file), and deactivate the option ```run_variant_calling```. $\textcolor{red}{\textsf{Be careful if you want to analyse a mix of bacterial species}}$, we recommend you to modify the ```reference_genome``` and deactivate the option ```run_variant_calling```.


## Running CamPype

1. Activate the CamPype's environment:
    ```bash
    $ conda activate campype
    ```
1. Go to CamPype's directory:
    ```bash
    $ cd your/path/to/CamPype
    ```
1. In case you want to run CamPype in the FASTQ mode, we encourage you to perform first a quality control step to check how good are your raw reads and adjust the read quality control filtering step (remember to include the path of the fastq files in the CamPype/input_files.csv file):
    ```bash
    $ campype_qc.py
    ```   
    A quality control analysis will be performed in each fastq file and a summary HTML report will be generated for fast interpretation inside the outputdirectoy of CamPype named fastq_quality_control. Check this [video](https://www.youtube.com/watch?v=bz93ReOv87Y) to know how to interpretate the results.
    
1. Once you have set the configuration, run CamPype:
    ```bash
    $ bash -i campype
    ```
1. \(*) You can deactivate the environment when you are finished:
    ```bash
    $ conda deactivate
    ```

## Output
The results of CamPype are stored in very detailed directories for each analysis, with separate subdirectories for each tool and isolate. The files will be generated for analysis tracking due to execution error. An interactive HTML summary report will be generated at the end of the analysis to simplify the task of data visualization and interpretation. This HTML file can be opened on any Web browser. An example of report can be found here.

You can generate the report after CamPype analysis by executing the following command in the Linux terminal:
* For raw fastq reads as input:
```
Rscript -e "rmarkdown::render('CamPype_Report_long.Rmd', params = list(directory = '~/path/to/data'))"
```
* For assembled genomes as input:
```
Rscript -e "rmarkdown::render('CamPype_Report_short.Rmd', params = list(directory = '~/path/to/data'))"
```
In both cases, you will have to change ```'~/path/to/data'``` with the corresponding path of the CamPype output directory containing the output files required to create the summary HTML report.


## How to update CamPype (Linux)

We recommend to update CamPype when newer versions are launched:

1. Make sure you are not in CamPype conda environment
    ```bash
    $ conda deactivate
    ```
1. Save a copy of your configuration files. Updating CamPype will overwrite you configuration files because this files properties may change with newer versions of CamPype.

1. Run ./updatecampype
    ```bash
    $ ./updatecampype
    ```
You should answer YES to the first question (Are you sure you want to remove your configuration files?) to update configuration files and YES to the second question (Proceed ([y]/n)?) to update all the packages and tools included in CamPype. This might take several minutes.


## How to uninstall CamPype (Linux)

1. Go to CamPype's directory
    ```bash
    $ cd your/path/to/CamPype
    ```
1. Make sure you are not in CamPype conda environment
    ```bash
    $ conda deactivate
    ```
1. The CamPype directory and every file within will be removed. Make sure you don't have anything important in this directory (results, data, etc).

1. Run ./uninstallcampype
    ```bash
    $ ./uninstallcampype
    ```

## FAQ
1. Prokka stops running with this error:
```
Could not run command: cat \/home\/CamPype_OUTPUT_20220511_131550\/Prokka_annotation\/NCTC11168\/NCTC11168\.IS\.tmp\.35844\.faa | parallel --gnu --plain -j 8 --block 313 --recstart '>' --pipe blastp -query - -db /home/instalador/anaconda3/envs/campype/db/kingdom/Bacteria/IS -evalue 1e-30 -qcov_hsp_perc 90 -num_threads 1 -num_descriptions 1 -num_alignments 1 -seg no > \/home\/CamPype_OUTPUT_20220511_131550\/Prokka_annotation\/NCTC11168\/NCTC11168\.IS\.tmp\.35844\.blast 2> /dev/null
```

Run ```prokka --setupdb``` first and execute CamPype again.

2.  ABRicate can't find any gen and this message appears: ```BLAST Database error: Error pre-fetching sequence data```

Run ```abricate --setupdb``` first and execute CamPype again.


## Citation
Please cite CamPype whenever you use it as:

Irene Ortega-Sanz, Jose A. Barbero and Antonio Canepa. CamPype. Available at [https://github.com/JoseBarbero/CamPype](https://github.com/JoseBarbero/CamPype)


## Contact:
For questions, bugs and suggestions, please open a new issue or contact us at iortega@ubu.es or jabarbero@ubu.es.
