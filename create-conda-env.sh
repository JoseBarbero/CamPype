export PATH=~/miniconda3/condabin:$PATH
conda config --append channels conda-forge
conda config --append channels bioconda
conda env create -f campype_env.yml
conda init
