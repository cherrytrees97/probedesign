# Probedesign 

## Table of Contents
[General Workflow](#general-workflow)

[Pre-Requisites](#pre-requisites)

[Preparing your alignment](#preparing-your-alignment)

A suite of Python scripts for designing TaqMan qPCR assays. 

## General workflow

TaqMan qPCR assays are designed in three steps: 

1. Assemble alignment and BLAST database (if checking for specificity)
2. Design and select viable TaqMan probe sequence 
3. Design and select viable primer sets for the chosen TaqMan probe

## Pre-Requisites

A conda environment `.yml` file is provided in this repo for easy setup. Otherwise, the scripts the following Python packages to be installed:
 
* `biopython`>=1.79
* `pandas`>=1.4.2
* `primer3-py`>=0.6.1

Additionally, `BLAST+`>=2.12 needs to be installed in your environment. 

## Preparing your alignment

The input for both the probe and primer design script is an alignment of the target loci sequences for your target species, in fasta format. 
```
>sequence_1
ATCGATC....ATGACATG
>sequence_2
ATCGATC....ATGACATG
>sequence_3
ATCGATC....ATGACATG
```

Generally, we use [MAFFT](https://mafft.cbrc.jp/alignment/software/) with the `--auto` flag to prepare our alignments. 

## Preparing a BLAST database for specificity checks

If you want to check the specificity of your probes and primers, you will need to create a custom BLAST database containing the sequences that you want to check.

Sequences for your BLAST database should be obtained from [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nucleotide/) in GenBank flat file format. Refer to the Appendix for more information on generating an Entrez query. Once you have downloaded the Genbank file, you will need to run two helper scripts: 
1. `process_genbank_db.py`: this script will output a FASTA file containing a sequence for each entry in the Genbank file, as well as parsing out relevant data into a metadata `.csv` file for easy reference. 
```
python ./process_genbank_db.py 
```
2. `generate_taxid_map.py`: this script will output a `.txt` file in tab delimited format that links the NCBI IDs to the tax IDs. 
```
python ./generate_taxid_map.py 
```

After running both of these scripts, you can use the following command to generate the database, using the FASTA and `.txt` file as inputs:

```
makeblastdb -in $fasta.fasta -parse_seqids -taxid_map $taxid_map.txt -dbtype nucl
```

Additionally, you will need to download the `taxdb.tar.gz` file provided by [NCBI](https://ftp.ncbi.nlm.nih.gov/blast/db/) in order for BLAST to recognize and convert tax IDs into species names. These files should be stored in the same directory as the database files. 

## Designing probes with probesearch.py
Ensure that you have activated the correct conda environment prior to running probesearch.py, and that you have the following files: 
1. FASTA alignment of all your target loci sequences
2. BLAST database (if checking specificity)
3. taxdb files (if checking specificity)

The program can be run with the following command (default settings) to find probes for the consensus sequence of the input alignment : 
```
python ./probesearch.py alignment_path --blastdb blastdb_path.fasta
```

If you do not want to perform specificity checking: 
```
python ./probesearch.py alignment_path --no_sens_spec
```

Probes can be filtered by two parameters: 
1. sequence representation (e.g. 50% of sequences in alignment contribute to the consensus sequence at the probe position)
2. annealing temperature

```
#Filter for sequence representation (above 70%) and annealing temperature (67 - 70 C)
python ./probesearch.py \
    alignment_path \
    --blastdb blastdb_path.fasta \
    --filter_seq_rep \
    --filter_min 0.6 \
    --filter_tm \
    --min_tm 67 \
    --max_tm 70
```

To increase the speed that the BLAST jobs are performed, you can pass a value equal to the number of processes you want to spawn with the flag `--mp_job`. 

```
#Spawn two processes when BLASTing probe sequences
python ./probesearch.py alignment_path --blastdb blastdb_path.fasta --mp_job 2
```

probedesign.py will output several files: 
1. `probe_candidates.csv` - contains the generated probes and their properties
2. `probe_blast/` - contains `.csv` files labeled by the probe root position and length that correspond to the BLAST results of that probe
3. `probesearch.log` - log file for the program

The probe that best fits your needs should be selected from the `probe_candidates.csv` file. Each probe has a predicted annealing temperature, and if enabled, a score between 0-2 that represents how sensitive and specific the probe is. The probe with the highest score and correct annealing temperature should be selected. 

## Designing primers using primersearch.py 
Ensure that you have activated the correct conda environment prior to running probesearch.py, and that you have the following files: 
1. FASTA alignment of all your target loci sequences
2. BLAST database (if checking specificity)
3. taxdb files (if checking specificity)

Additionally, `probesearch.py` should have been run and you should have selected a couple probes that you are interested in designing primers for. Note the putative probes' root positions and length. 

The program can be run with the following command (default settings) to find primers for the consensus sequence of the input alignment for the input probe parameters: 

```
python ./primersearch.py alignment_path probe_root probe_len --blastdb blastdb_path.fasta
```

If you do not wish to check specificity: 
```
python ./primersearch.py alignment_path probe_root probe_len --no_sens_spec
```

Primers can be filtered by three parameters: 
1. sequence representation  (e.g. 50% of sequences in alignment contribute to the consensus sequence at the primer positions) (default: 0.5)
2. annealing temperature (default: 55.0 to 63.0)
3. max annealing temperature difference between the forward and reverse primer (default: 5.0)

```
#Filter for annealing temperature between 55.0 to 63.0, with a sequence representation over 50% and a max annealing temperature difference of 5.0
python ./primersearch.py \
    alignment_path \
    probe_root \
    probe_len \
    --blastdb blastdb_path.fasta \
    --filter_tm \
    --max_primer_tm 55 \
    --min_primer_tm 63 \
    --max_tm_diff 5 \
    --filter_seq_rep \
    --filter_min 0.6 \
```

To increase the speed that the BLAST jobs are performed, you can pass a value equal to the number of processes you want to spawn with the flag `--mp_job`. 

```
python ./primersearch.py alignment_path probe_root probe_len --blastdb blastdb_path.fasta --mp_job 2
```

primerdesign.py will output several files: 
1. `primer_pairs.csv` - contains the generated pairs and their properties
2. `fw_blast/` - contains `.csv` files labeled by the forward primer root position and length that correspond to the BLAST results of that primer
3. `rev_blast/` - contains `.csv` files labeled by the reverse primer root position and length that correspond to the BLAST results of that primer
3. `primersearch.log` - log file for the program

The primer pair that best fits your needs should be selected from the `probe_candidates.csv` file. Primer pairs will have a predicted annealing temperature calculated using primer3 with default settings for each primer. There are individual scores for each primer, as well as a combined score for the primer pair, which takes into account only accessions that are found in the BLAST results of both primers. The optimal primer pair will have a predicted annealing temperature roughly 10 degrees less than the probe annealing temperature or around 60 C, and will have a high combined score. 

## Appendices
### Appendix I: Generating an Entrez Query
See 'general-scripts' wiki [entry](Nothing here yet). 
### Appendex II: probesearch.py arguments
```
probesearch.py - identify viable probes in an alignment for given target sequences

positional arguments:
  target_alignment_path
                        Path to target alignment file, fasta format

optional arguments:
  -h, --help            show this help message and exit

Probe parameters:
  --target_start target_start
                        Start coordinate of target region, 1-based coordinates
  --target_end target_end
                        End coordinate of target region, 1-based coordinates
  --min_probe_len MIN_PROBE_LEN
                        Minimum probe length (default=17)
  --max_probe_len MAX_PROBE_LEN
                        Maximum probe length (default=22)

BLAST parameters:
  --no_sens_spec_check  Flag to not check the putative probes for their specificity and sensitivity
  --blastdb BLASTDB     Name of blastdb
  --mp_job NUM_JOBS, -m NUM_JOBS
                        Number of processes to spawn to handle BLAST jobs. (Default=1)

Output parameters:
  --output output_directory, -o output_directory
                        Output path

Filter parameters:
  --filter_seq_rep, -fs
                        Filter by probes returned by sequence representation
  --filter_min MIN_SEQ_REP
                        Minimum percentage of sequences that need to be represented for probes to be returned. Default = 0.5
  --filter_tm, -ft      Filter probes by tm range
  --min_tm MIN_TM       Minimum tm of probe for it to be returned. Default = 68.0
  --max_tm MAX_TM       Maximum tm of probe for it to be returned. Default = 70.0
```

### Appendix III: primersearch.py arguments
```
primersearch.py - identify viable primers for a given probe

positional arguments:
  target_alignment_path
                        Path to target alignment file, fasta format
  pb_start              Start coordinate of probe, 1-based coordinates
  pb_len                Length of the probe

optional arguments:
  -h, --help            show this help message and exit

Primer parameters:
  --min_primer_len MIN_PRIMER_LEN
                        Minimum primer length
  --max_primer_len MAX_PRIMER_LEN
                        Maximum primer length

BLAST parameters:
  --no_sens_spec        Flag to not check the putative probes for their specificity and sensitivity
  --blastdb BLASTDB     Name of blastdb
  --mp_job NUM_JOBS, -m NUM_JOBS
                        Number of processes to spawn to handle BLAST jobs. (Default=1)

Output parmaeters:
  --output_path OUTPUT_PATH, -o OUTPUT_PATH
                        Output path

Output Filter parameters:
  --filter_tm, -ft      Filter primers by tm diff, as well as min and max tm
  --max_tm_diff-d maximum_tm_diff
                        Maximum temperature difference between forward and reverse
  --max_primer_tm MAX_TM
                        Maximum tm for primer.
  --min_primer_tm MIN_TM
                        Minimum tm for primer.
  --filter_seq_rep, -fs
                        Filter by probes returned by sequence representation
  --filter_min MIN_SEQ_REP
                        Minimum percentage of sequences that need to be represented for probes to be returned. Default = 0.5
```
