# prbh

Python program for detecting Phospho-Regulated Basic Hydrophobic motifs
in proteins as described in [Bailey and Prehoda](https://doi.org/10.1016/j.devcel.2015.09.016)

## Usage

On the command line give the FASTA file containing protein
sequences (can be gzipped) and an optional output file (otherwise output is sent to standard output).

```
./prbh.py fasta_file.faa -o output.txt
```
