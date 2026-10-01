[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23046565.svg)](https://doi.org/10.5281/zenodo.23046565)
# Coming soon!
MotifAI is still under active development, but will come online soon.
# MotifAI: A machine learning framework for prioritizing transcription factor binding motifs in putative regulatory DNA.
Current software, such as the popular motif scanning tool Find Individual Motif Occurrences (FIMO), can identify all possible transcription factor (TF) binding motifs in input DNA sequences of interest. However, even from relatively short input sequences, these tools regularly output overwhelming quantities of motifs, making manual selection of high-confidence motif "hits" cumbersome and uncertain work. Trained on a corpus of genetic, epigenetic, and evolutionary datasets in the model plant species *Arabidopsis thaliana*, MotifAI is an ML-powered tool for ranking/prioritizing identified motifs based on likelihood of functionality in regulating cognate gene expression.
# Installation & Configuration
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
## Configuration
Configure your environment for MotifAI with the configuration script. Run:
```bash
motifai_Configure
```
This requires a recent version of mamba or conda to be installed and set on your system path. It creates three conda environments necessary for downloading data and running MotifAI.
# Data Download
The model underlying MotifAI is an XGBoost classification model trained on a feature table consisting of 28 total features. As a result, a feature table must be constructed based on the user's input DNA sequences/coordinates, which the trained model weights uses to. Thus, the *Arabidopsis thaliana* data files (e.g., CNS coordinates, genome-wide DNA methylation coverages, etc.) used to construct the original feature table must be downloaded from a Zenodo repository, and the trained model weights downloaded from a HuggingFace repository.

To download the necessary data & model weights, simply run:
```bash
motifai_FetchData -o /path/to/desired/output_directory/
```
This outputs in the specified directory:
1. `arabidopsis-motifai-v1.0-xgbc_FeatureTable.rds`: The original feature table used to train the model underlying MotifAI, in RDS format.
2. `Zenodo_files/`: Directory holding all data files necessary for feature table construction.
3. `model/`: Directory holding trained model weights.

Now you are ready to use MotifAI!
# Usage
Command-line usage of MotifAI is simple but requires three arguments:
1. `-i`: DNA sequences in which you desire to identify and rank motifs.
2. `-z`: Path to the `Zenodo_files` (determined in the last step)
3. `-o`: Path you want MotifAI's output to save to.
```bash
motifai -i dna_sequences.fa \
		-z /path/to/Zenodo_files/
		-o /path/to/output_directory/
```
