# CamPype

CamPype is a pipeline for the analysis of Illumina paired-end sequencing data and/or whole bacterial genomes. The development of the workflow is mainly intended for the analysis of <em>Campylobacter jejuni/coli</em> genomes, although any other bacterial genus can be analyzed as well. CamPype is specially designed for users without knowledge of bioinformatics or programming, so the ease of installation and execution are the fundamentals of its development. Moreover, CamPype is a user-customizable workflow that allows you to select the analysis and the tools you are interested in. 


Here you wil find the schema of CamPype. CamPype allows the user to previously check the quality of sequencing raw data in an independent step to optimize the read filtering analysis. Moreover, bacteria identification can be performed on the filtered fastq reads when raw reads are provided or after genome assembly when contigs are used.

![scheme](https://github.com/JoseBarbero/CamPype/assets/58036036/1589d4a8-cd1a-42eb-be38-f2256086272f)
Software or databases are indicated in boxes, while discontinuous boxes indicate tools that users can deactivate.


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

1. After conda environments were created, update the databases of AMRFinder, Prokka and ABRicate:
    ```
    $ conda activate campype
    $ amrfinder -u
    $ prokka --setupdb
    $ abricate --setupdb
    ```
1. Additionally, CamPype allows you to check for read contamination and determine bacteria taxonomy using Kraken2. The installation of Kraken2 is optional to avoid possible storage limitations as it requires the use of a heavy database that requires high free disk space, but it is not needed for Campype if you are not interested in this analysis. If you want to install this module, at least 8 GB will be occupied to store the MiniKraken_8GB_202003 database. This database is enough for bacteria identification and CamPype performance, but if you have enough disk space, you can download and install "Standard Kraken2 Databases" following the instructions [here](https://lomanlab.github.io/mockcommunity/mc_databases.html) for better sensitivity. To install Kraken2 in CamPype run:
    ```bash
    $ install_kraken.sh
    ```

## Set input files and configuration

CamPype can run on two modes depending on the input files. The FASTQ mode analyses (un)compressed raw reads in fastq format, while the FASTA mode analyses assembled genomes in fasta format.

1.  Before running CamPype, you must indicate the location of the input files in the CamPype/input_files.csv file (TAB as separator). For fastq files, you need to indicate the path of each pair of reads in the Forward and Reverse columns, while for fasta files both the Forward and Reverse columns will refer to the path of the assembled genome, that is, the content of both columns will be the same. We recommend you to indicate the Genus and Species of your genomes in case you known this information beforehand (this taxonomy will be considered even in you run the bacteria identification module). Please, indicate the species in case you indicate the genus. If not, we strongly encourage you to activate the bacteria identification analysis (as below explained).

    | Samples        | Forward           | Reverse  | Genus  | Species  |
    | ------------- |:-------------:|:-----:|:-----:|:-----:|
    | sample_ID  | /path/to/your/forward/fastq1_file.fastq | /path/to/your/reverse/fastq1_file.fastq | YourStrainGenus | YourStrainSpecies
    | sample_ID  | /path/to/your/forward/fastq2_file.fastq | /path/to/your/reverse/fastq2_file.fastq | YourStrainGenus | YourStrainSpecies
    | sample_ID  | /path/to/your/forward/fastq3_file.fastq | /path/to/your/reverse/fastq3_file.fastq | YourStrainGenus | YourStrainSpecies

    This structure must be respected anyway in that file. Make sure headers haven't changed and samples ID do not contain the dot symbol ```.``` $\textcolor{red}{\textsf{Be careful with typos!!!}}$
    
1. Set the modules you want to run in the CamPype/workflow_config.py file. There you will set your own running parameters for each tool and the select your tools of interest when possible.

1. Default settings are configured for <em>Campylobacter jejuni/coli</em>. If you want to use a different bacteria, we strongly recommend you to adapt the configuration of CamPype as previously explained. In particular, you must modify the ```reference_genome```, deactivate the option ```include_cc```, and use abricate for virulence genes searching or/and use your own virulence genes database with BLAST instead (indicate this accordingly in the CamPype/workflow_config.py file). $\textcolor{red}{\textsf{Be careful if you want to analyse a mix of bacterial species}}$, we recommend you to delete a ```reference_genome``` and deactivate the option ```run_variant_calling```.


## Test CamPype
An optional test can be run to check the correct installation of CamPype.
1. Activate the CamPype's environment:
    ```bash
    $ conda activate campype
    ```
1. Go to CamPype's directory:
    ```bash
    $ cd your/path/to/CamPype
    ```
1. Run the CamPype's test:
    ```bash
    $ campype_test.sh
    ```   
An example of the HTML report you will get can be found here.


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
    $ bash -i campype_qc
    ```   
    A quality control analysis will be performed in each fastq file and a summary HTML report will be generated for fast visualization in the directory fastq_quality_control, that will be located inside the output directoy of CamPype named as you indicated in the CamPype/workflow_config.py file. Check this [video](https://www.youtube.com/watch?v=bz93ReOv87Y) to know how to understand these results.
    
1. Once you have set the configuration, run CamPype:
    ```bash
    $ bash -i campype
    ```
    The results will be located in the output directoy of CamPype named as you indicated in the CamPype/workflow_config.py file. If you want to store these results in the same directory where the quality control analysis data are, remember to indicate that directory in the configuration file. 
    
1. You can deactivate the environment when you are finished:
    ```bash
    $ conda deactivate
    ```

## Output
The results of CamPype are stored in very detailed directories for each analysis, with separate subdirectories for each tool and isolate. The files will be generated for analysis tracking due to execution error. An interactive HTML summary report will be generated at the end of the analysis to simplify the task of data visualization and interpretation. This HTML file can be opened on any Web browser. Examples of reports can be found here:
* [Analysis with 5 Campylobacter jejuni and 5 Campylobacter coli (raw reads)](https://josebarbero.github.io/CamPype/example_report/CamPype_Report_long_first_case_study.html)
* [Analysis with 44 Escherichia coli (assembled genomes)](link)

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
$\textcolor{red}{\textsf{Be careful, the CamPype directory and every file within will be removed.}}$ 
$\textcolor{red}{\textsf{Make sure you don't have anything important in this directory (results, data, etc).}}$

1. Run ./uninstallcampype
    ```bash
    $ ./uninstallcampype
    ```

## FAQ
* Prokka stops running with this error:
  ```
  Could not run command: cat \/home\/CamPype_OUTPUT_20220511_131550\/Prokka_annotation\/NCTC11168\/NCTC11168\.IS\.tmp\.35844\.faa | parallel --gnu --plain -j 8 --block 313 --recstart '>' --pipe blastp -query - -db /home/instalador/anaconda3/envs/campype/db/kingdom/Bacteria/IS -evalue 1e-30 -qcov_hsp_perc 90 -num_threads 1 -num_descriptions 1 -num_alignments 1 -seg no > \/home\/CamPype_OUTPUT_20220511_131550\/Prokka_annotation\/NCTC11168\/NCTC11168\.IS\.tmp\.35844\.blast 2> /dev/null
  ```
  Activate the CamPype's directory ```conda activate campype```, run ```prokka --setupdb``` first, and execute CamPype again.

* ABRicate can't find any gen and this message appears: ```BLAST Database error: Error pre-fetching sequence data```

  Activate the CamPype's directory ```conda activate campype```, run ```abricate --setupdb``` first, and execute CamPype again.

* If you find missing font problems running Mauve, you should install the required fonts:
  ```bash
  $ sudo apt-get install ttf-dejavu
  ```


## Citation
Please cite CamPype whenever you use it as:

Irene Ortega-Sanz, Jose A. Barbero and Antonio Canepa. CamPype. Available at [https://github.com/JoseBarbero/CamPype](https://github.com/JoseBarbero/CamPype)


## Contact
For questions, bugs and suggestions, please open a new issue or contact us at iortega@ubu.es or jabarbero@ubu.es.
