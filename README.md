ScFvFinder is a tool for identifying and finding valid ScFv sequences from long read sequencing from the fastq file.

The general pipeline and tools used:

NanoFilt
Cutadapt
awk
blastn - needs to be installed by user
igblast - if not installed then will be installed into the data_files directory


The tool assumes that you have variable chains separated by a linker, you supply a reference fasta file that includes the linker DNA sequence, e.g.:

'''
>Linker
GGCGGCGGCGGCAGC
'''

blastn is then used to identify sequences which contain the linker and then contain sequences pre and post the linker that fit with the likely size of a complete ScFv

these "pre" and "post" linker sequences are then separated.

igblast is then used to identify the regions of the variable chains.

This data is then processed to identify full variable chain sequences that contain fwr1, cdr1, fwr2, cdr2, fwr3, cdr3 and fwr4

Stop codons within the fwr regions are then mutated out to the germline reference as they are most likely sequencing artefacts.

Sequences with stop codons in the cdrs are excluded

The sequences are then put back together, orientation calculated (H-L or L-H) and then a final check to ensure no stop codons have gone through.

