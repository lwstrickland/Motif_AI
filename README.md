[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23046565.svg)](https://doi.org/10.5281/zenodo.23046565)
# MotifAI: A machine learning framework for prioritizing transcription factor binding motifs in putative regulatory DNA.
Current software, such as the popular motif scanning tool Find Individual Motif Occurrences (FIMO), can identify all possible transcription factor (TF) binding motifs in input DNA sequences of interest. However, even from relatively short input sequences, these tools regularly output overwhelming quantities of motifs, making manual selection of high-confidence motif "hits" cumbersome and uncertain work. Trained on a corpus of genetic, epigenetic, and evolutionary datasets in the model plant species *Arabidopsis thaliana*, MotifAI is an ML-powered tool for ranking/prioritizing identified motifs based on likelihood of genuine TF-binding potential.
# Installation & Configuration
## Installation (Linux)
Download MotifAI via `git clone`:
```bash
git clone https://github.com/lwstrickland/Motif_AI.git
```
This will create a new directory in your working directory called `Motif_AI`.
Then, simply add the executable scripts to your user path, either by running the following on the command line (required once per session) **OR** by adding the following to your `.bashrc` or `.zshrc` (one-time only):
```bash
export PATH="$PATH:/your/path/to/Motif_AI/scripts"
```
To test if MotifAI is now accessible from wherever on the command line, type in the commands to bring up their respective help pages:
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
**NOTE:** This step requires an Internet connection.

The model underlying MotifAI is an XGBoost classification model trained on 28 total features. As a result, a feature table must be constructed based on the user's input DNA sequences/coordinates, which the model uses to make predictions on identified motifs. Thus, the *Arabidopsis thaliana* data files (e.g., CNS coordinates, genome-wide DNA methylation coverages, etc.) used to construct the feature table must be downloaded from a Zenodo repository, and the trained model weights downloaded from a HuggingFace repository.

To download the necessary data & model, simply run:
```bash
motifai_FetchData -o /path/to/data/
```
This outputs in the specified directory `/path/to/data/`:
1. `arabidopsis-motifai-v1.0-xgbc_FeatureTable.rds`: The original feature table used to train the model underlying MotifAI, in RDS format.
2. `Zenodo_files/`: Directory holding all data files necessary for feature table construction.
3. `model/`: Directory holding trained model weights.

Now you are ready to use MotifAI!
# Usage
Command-line usage of MotifAI is simple but requires three arguments:
1. `-i`: DNA sequences in which you desire to identify and rank motifs.
2. `-d`: Path to directory holding `Zenodo_files/` and `model/` subdirectories (defined by you by `motifai_FetchData -o /path/to/data/`).
3. `-o`: Path you want MotifAI's output to save to.
```bash
motifai -i dna_sequences.fa \
		-d /path/to/data/
		-o /path/to/output_directory/
```
The primary output of MotifAI is `motifai_Results_mm.dd.yyyy.tsv`:
```
motif	sequence	TF	probability	rank
MA2007.2_Chr3:433-440	CACCAAAC	MYB107	0.81969	1
MA1675.2_Chr3:12770495-12770502	ACGCAACT	NAC029	0.69437	2
MA1245.1_Chr3:552-566	ACAGCAGCACCGTAG	ERF112	0.6531	3
UN0846.1_Chr3:433-441	GTTTGGTGA	AT4G26030	0.62249	4
```
1. `motif`: JASPAR motif code + genomic coordinates of motif
2. `sequence`: Motif sequence
3. `TF`: Predicted motif's corresponding transcription factor
4. `probability`: Probability of motif's genuine TF-binding potential, assigned by MotifAI
5. `rank`: Motif's rank/priority

Other outputs includes:
- `FIMO_results/`: Raw output of running FIMO on input DNA sequences
- `feature_table.rds` and `feature_table.tsv`: Constructed feature table for identified motifs in RDS and TSV formats
## Notes
**1)** In order for MotifAI to properly parse the genomic coordinates of your input DNA sequences (`dna_sequences.fa`), the FASTA headers must contain the genomic coordinates, like this:
```bash
>Chr1:456765-456900
CAGATCATTTA . . .
```
This is the default output format for `bedtools getfasta`, which extracts coordinate-defined DNA sequences from FASTA files.

**2)** Currently, MotifAI is only built to rank DNA sequences in *Arabidopsis thaliana* (TAIR10 genome build). The developers are currently working on enabling broader functionality for motif ranking in other important plant species. Stay tuned!
# Contact
The MotifAI developers welcome suggestions for making the tool work better for the gene regulation community. If you wish to offer such a suggestion, please do not hesitate to reach the developers directly: stric132@msu.edu. Thank you for using MotifAI!