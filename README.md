# ScFvFinder

ScFvFinder is a tool for identifying and finding valid ScFv sequences from long read sequencing from the fastq file.


To use:

## Installation

For now only MacOS is supported, linux support is imminent, windows support unknown

These commands only need to be run once

clone the repo into an appropriate directory

```bash
cd ~/Tools # Just for example
git clone https://github.com/bend91/ScFvFinder.git
cd ScFvFinder
```
Optional but recommended - make a virtual enviornment, recommended to use python version 3.10 or higher, there are some broken dependencies in cutadapt with python 3.7, not sure where they are fixed

To check python versions available
```bash
whereis python3
```
```bash
python3.10 -m venv ~/VENV/ScFvFinder
source ~/VENV/ScFvFinder/bin/activate
```
Install the python requirements (note this doesn't install igblast or blastn)
```bash
pip install -r requirements.txt
```

## Using the tool
```bash
python scfv_find.py
```

---

The general pipeline and tools used:

- NanoFilt (https://github.com/wdecoster/nanofilt)
- Cutadapt (https://github.com/marcelm/cutadapt)
- blastn - needs to be installed by user
- igblast - if not installed then will be installed into the data_files directory





The tool assumes that you have variable chains separated by a linker, you supply a reference fasta file that includes the linker DNA sequence, e.g.:

```fasta
>Linker
GGCGGCGGCGGCAGC
```


blastn is then used to identify sequences which contain the linker and then contain sequences pre and post the linker that fit with the likely size of a complete ScFv

these "pre" and "post" linker sequences are then separated.

igblast is then used to identify the regions of the variable chains.

This data is then processed to identify full variable chain sequences that contain fwr1, cdr1, fwr2, cdr2, fwr3, cdr3 and fwr4

Stop codons within the fwr regions are then mutated out to the germline reference as they are most likely sequencing artefacts.

Sequences with stop codons in the cdrs are excluded

The sequences are then put back together, orientation calculated (H-L or L-H) and then a final check to ensure no stop codons have gone through.

