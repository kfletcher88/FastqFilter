A wrapper script to filter fastq reads against a custom nucleotide database and extract:
	1. Reads which do not map to anything
	2. Reads which map to reference organisms of interest.

To use this script, the user must provide a database file where the target reference sequences are easily identifiable by a common string.

To do this, users may want to rename the sequesnces in target organism fasta files before collating them in to a larger file with non-target organisms.

A simple `awk` one liner will do the trick:

```
awk -v '/^>/{print ">"TargetSeq"" ++i; next}{print}' < Input.fasta > Output.fasta
```

