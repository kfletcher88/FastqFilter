A wrapper script to filter fastq reads against a custom nucleotide database and extract:

1. Reads which do not map to anything

2. Reads which map to reference organisms of interest.

These read files will be placed in a new directory, named by the user with option `-o`

The software required by this wrapper are:

1. [BWA](https://arxiv.org/abs/1303.3997v2 "Li et al. 2013")
2. [SAMtools](https://www.ncbi.nlm.nih.gov/pubmed/19505943 "Li et al. 2009")
3. [PairFQ_lite](https://github.com/sestaton/Pairfq "Staton GitHub repo")

Upon initiation the script will check that these executables are in the PATH.

For usage type:

    ./FastqFilter.sh -h

To use this script, the user must provide a database file where the target reference sequences are easily identifiable by a common string.

To do this, users may want to rename the sequesnces in target organism fasta files before collating them in to a larger file with non-target organisms.

A simple `awk` one liner will do the trick:

```
awk -v '/^>/{print ">""TargetSeq-" ++i; next}{print}' < Input.fasta > Output.fasta
```

This will output a new fasta file, where every sequence is named "TargetSeq-[num]" starting at 1.


The shell script wraps four basic commands and adds in some safe-guards to ensure that nothing important gets accidently deleted.

The basic commands ran are:

```
bwa index [Database.fa] 
bwa mem [Database.fa] [Read1.fq] [Read2.fq] | samtools view -bT [Database.fa] - | samtools sort -o [Prefix].sorted.bam -

samtools index [Prefix].sorted.bam
samtools fastq -f 4 -1 [Prefix].UnMapped.1.fq -2 [Prefix].UnMapped.2.fq [Prefix].sorted.bam
samtools view [Prefix].sorted.bam | awk -v var="[Identifier]" '$3 ~ var' | samtools view -bT [Database.fa] | samtools fastq -1 [Prefix].Mapped.1.fq -2 [Prefix].Mapped.2.fq -

pairfq_lite.pl addinfo -i [Prefix].Mapped.1.fq -o [Prefix].Mapped.info.1.fq -p 1
pairfq_lite.pl addinfo -i [Prefix].Mapped.2.fq -o [Prefix].Mapped.info.2.fq -p 2
pairfq_lite.pl addinfo -i [Prefix].UnMapped.1.fq -o [Prefix].UnMapped.info.1.fq -p 1
pairfq_lite.pl addinfo -i [Prefix].UnMapped.2.fq -o [Prefix].UnMapped.info.2.fq -p 2

pairfq_lite.pl makepairs -f [Prefix].Mapped.info.1.fq -r [Prefix].Mapped.info.2.fq -fp [Prefix].Mapped.RePair.1.fq -rp [Prefix].Mapped.RePair.2.fq -fs [Prefix].Mapped.NoPair.1.fq -rs [Prefix].Mapped.NoPair.2.fq
pairfq_lite.pl makepairs -f [Prefix].UnMapped.info.1.fq -r [Prefix].UnMapped.info.2.fq -fp [Prefix].UnMapped.RePair.1.fq -rp [Prefix].UnMapped.RePair.2.fq -fs [Prefix].UnMapped.NoPair.1.fq -rs [Prefix].UnMapped.NoPair.2.fq
```

Where the text contained in square brackets are variables provided to the shell script. Specifically:

`-d` Database.fa

`-1` Read1.fq

`-2` Read2.fq

`-o` Prefix

`-s` Identifier string, by which sequences, closely related to the organism under study, can be identified for positive filtering.

In the end, four fastq files are output and are in their correct pairs.

Running this shell script will also automatically clean up intermediate files. If the commands are ran manually, as above, then clean up is not automatic. 


