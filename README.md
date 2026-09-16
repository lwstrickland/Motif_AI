# Coming soon!
MotifAI is still under active development, but will come online soon.
# MotifAI: A machine learning framework for prioritizing transcription factor binding motifs in putative regulatory DNA.
Current software, such as the popular motif scanning tool Find Individual Motifs Occurrences (FIMO), is able to identify all possible transcription factor (TF) binding motifs in input DNA sequences of interest. However, even from relatively short input sequences, these tools regularly output overwhelming quantities of motifs, making manual selection of high-confidence motif "hits" cumbersome and uncertain work. Trained on a corpus of genetic & epigenetic data in the model plant species *Arabidopsis thaliana*, MotifAI is an ML-powered tool for ranking/prioritizing identified motifs based on likelihood of functionality in regulating cognate gene expression.
# Installation & Setup
## Installation
Download MotifAI via `git clone`:
```bash
git clone git@github.com:lwstrickland/Motif_AI.git
```
This will create a new directory in your working directory called `Motif_AI`.
Then, simply add executable scripts to your user path, either by running the following line on the command line (required once per session) by adding the following lines to your `.bashrc` or `.zshrc` (one-time only):
```bash
export PATH="$PATH:/your/path/to/Motif_AI"
```
To test if MotifAI functionality is now accessible from wherever on the command line, type in the commands to bring up their respective help pages:
```bash
motifai_FetchData # executable for fetching necessary data
motifai # executable for running MotifAI
```
## Setup
For simplicity, it is recommended to setup the required environment and install the necessary softwares for MotifAI using conda. Create ***two conda environments***: one for identifying motifs of present in your input DNA sequences using FIMO, a tool distributed as a part of the MEME suite, and another for prioritizing those motifs with MotifAI.
	This is necessary because MEME & everything else do not play nicely in the same conda environment (i.e., the environment simply refuses to solve due to package version incompatibilities).
#### Conda env: `meme`
Create a conda env for running FIMO; this can be done with two simple commands:
```bash
conda create -n meme -c conda-forge -c bioconda meme=5.5.4
conda activate meme
fimo # test
```
If this pulls up the FIMO help page on the screen, you should be good to go.
#### Conda env: `motifai`
The MotifAI software is implemented using helper scripts in the R programming language. If you do not R installed, you can simply use this conda env setup through a `environment.yml` file, part of the official GitHub distribution for MotifAI. This file is downloaded upon `git clone` (contents shown below):
```yaml
name: motifai
channels:
  - conda-forge
  - bioconda
  - nodefaults
dependencies:
  - r-base=4.3
  - r-tidyverse
  - r-cowplot
  - r-showtext
  - r-showtextdb
  - r-sysfonts
  - r-caret
  - r-prroc
  - r-argparse
  - r-xgboost=1.7.6
  - bioconductor-genomicranges=1.54.1
  - bioconductor-rtracklayer=1.62.0
  - zenodo_get
  - huggingface_hub
```

Then, create the environment with the necessary packages and such:
```bash
conda env create -f /path/to/Motif_AI/environment.yaml
conda activate motifai
```
To save a written record of all packages & package versions installed in your conda environment:
```bash
conda env export --no-builds > environment-lock.yaml
```
