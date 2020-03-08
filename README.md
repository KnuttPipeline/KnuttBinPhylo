# KnuttBinPhylo

The bin (metagenome-assembled genomes/MAG) phylogenetic analysis pipeline

This is the version of the pipeline as it was used for the paper. Some data outputs weren't used and have been disabled.

## Setup

This pipeline requires a lot of disk space for its reference databases (~200GB)

1. Install the latest version of a conda distribution and Snakemake (>=5.10)
2. Clone this repository (paper branch)
3. Copy the bins into the `input` folder (ending in `.fa`)
4. Run the `paper` rule with conda enabled

``` sh
conda create -y -n snake -c bioconda -c conda-forge snakemake>=5.10
conda activate snake
git clone -b paper https://github.com/KnuttPipeline/KnuttBinPhylo.git
cd KnuttBinPhylo

snakemake -prj 16 --use-conda paper 2>&1 | tee run.log
```
