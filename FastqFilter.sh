#!/bin/bash
#Wrapper to map paired end next generation sequencing reads to a custom database and filter reads.

#MIT liscense
#Copyright 2017 Kyle Fletcher
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction,
#including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
#subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
#IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
#OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

usage="($basename "$0") [-h] [-d] [-1] [-2] [-M] [-o] [-s] [-T] [-L] --  A wrapper to run align paired end next generation sequencing reads to a custom database and output:
	1. A pair of files containing paired, unmapped reads
	2. A pair of files containing paired reads which mapped to target reference sequences

	These four files will be placed in a directory named with option o

        Please note, the database must be formtatted such that each target reference sequence is labelled with a custom string, to allow the identification of reads mapping to these sequences.
	The README.md contains instructions for how to do this, if required:

	Before combining all reference sequences together
        We have hard-coded in the blast word size (16) and number of threads (1). If you wish to change these then please edit the command on line 156

Options:
        -h show this help message
        -d bwa indexed fasta file - May be gzipped
        -1 Paired end read file 1 - May be gzipped
        -2 Paired end read file 2 - May be gzipped
        -M Single ended reads - May be gzipped
	-o Output prefix
        -s Unique string, by which target reference sequences may be identified
	-T Number of threads (default = 8)
	-L Dependency list
"

DepL="Dependencies required to run this script successfully are:

	Samtools v1.3 +			[http://samtools.sourceforge.net/]
	bwa mem (tested with v0.7.12) 	[https://sourceforge.net/projects/bio-bwa/files/]
	Pairfqlite.pl			[https://github.com/sestaton/Pairfq]

All executable should be located in the PATH
"
#How about an option to obtain only unmapped reads
#How about an option to obtain unpaired reads
while getopts ':hd:1:2:M:o:s:T:L' option; do
        case "$option" in
                h)  echo "$usage"
                         exit
                         ;;
                d)  DB=$OPTARG
                         ;;
                1)  R1=$OPTARG
                         ;;
                2)  R2=$OPTARG
                         ;;
		M)  SE=$OPTARG
			 ;;
                o)  Prefix=$OPTARG
                         ;;
                s)  string=$OPTARG
                         ;;
		L)  echo "$DepL"
			 exit
			 ;;
		T)  Threads=$OPTARG
			 ;;
                \?) printf "illegal option: -%s\n\n" "$OPTARG" >&2
                    echo "$usage"
                    exit 1
                         ;;
        esac
done


#Check Arguments provided
if [[ -z $DB || -z $Prefix || -z $string ]]; then
       echo "
       ERROR: Some required variables are not set
"
       echo "$usage"
       exit 1
fi

if [[ -z $R1 || -z $R2 || -z $SE ]]; then
	echo "
	ERROR: No read files provided
	Please provide either:
	Paired End files with flags -1 & -2
	Single End files with flag -M
"
fi

if [[ -n $SE ]]; then 
	echo "
	Single end reads provided
"
fi

if [[ -n $R1 && -n $R2 ]]; then
	echo "
	Paired end reads provided
"
fi

#Set defaults if not arguments provided
if [[ -z $Threads ]]; then
Threads=8
echo "
Threads not specified, will proceed with default [8]
"
fi

#Create Output directory
mkdir -p $Prefix-RFA
echo "Results will be output to directory:"
echo $PWD/$Prefix-RFA 

#Check for previously generated files to make sure they are not overwrote

#Check dependencies
SAM=$(samtools 2> SAMtmp ; grep 'Version:' SAMtmp | awk '$2 >= 1.3 {print "Samtools version", $2, "detected, ok to proceed"} ; $2 < 1.3 {print "Samtools version", $2, "detected. We required Samtools version >= 1.3"}')
#clean up
rm SAMtmp
if [[ $SAM =~ 1.3$ ]]; then
        echo "$SAM"
        exit 1
else
        echo "$SAM"
fi

BWArun=$(command -v bwa mem >//dev/null 2>&1 || { echo >&2 "bwa mem not detected, please check it is in your PATH"; })
if [[ $BWArun != "" ]]; then
        echo "$BWArun"
        exit 1
else
        echo "bwa mem detected"
fi
PFQ=$(command -v pairfq_lite.pl >//dev/null 2>&1 || { echo >&2 "Pairfq_lite.pl not detected. Please check that it is in your PATH"; })
if [[ $PFQ != "" ]]; then
	echo "$PFQ"
	exit 1
else
	echo "pairfq_lite.pl detected"
fi

#Check number of target sequences in the file
if [ -e $DB ]; then
echo "Counting number of target sequences present in the database"
	if [[ $DB =~ .gz$ ]]; then
	RefCount=$(zcat $DB | tee $DB.unzip | grep '^>' | grep -c $string -)
	else
	RefCount=$(grep '^>' $DB | grep -c $string -)
	fi
	echo "We detected $RefCount"
	if [[ $RefCount == 0 ]]; then
	echo "Reference sequences are required for positive mapping. Please check that you have formatted the database correctly"
	exit 1
	else
	echo "Moving on to Mapping"
	fi
else
echo "The database fasta file does not seem to be present, we will require this file for working with Samtools"
exit 1
fi

if [[ -e $DB.amb || -e $DB.ann || -e $DB.bwt || -e $DB.pac || -e $DB.sa ]] ; then
echo "All index files detected"
else
echo "Index files not detected, let's build them!"
	if [ -e $DB ] ; then
	bwa index $DB
	else
	echo "Database fasta file not present... we can't proceed. Please make sure the file is where you say it is"
	exit 1
	fi
fi

if [[ $DB =~ .gz$ ]]; then
	if [[ -n $R1 && -n $R2 ]]; then
	bwa mem -t $Threads $DB $R1 $R2 | samtools view -bT $DB.unzip - | samtools sort -o $Prefix.sorted.bam -
	samtools index $Prefix.sorted.bam
	samtools fastq -f 4 -1 $Prefix.UnMapped.1.fq -2 $Prefix.UnMapped.2.fq $Prefix.sorted.bam
	samtools view $Prefix.sorted.bam | awk -v var="$string" '$3 ~ var' | samtools view -bT $DB.unzip | samtools fastq -1 $Prefix.Mapped.1.fq -2 $Prefix.Mapped.2.fq -
	fi
	if [[ -n $SE ]]; then
	bwa mem -t $Threads $DB $SE | samtools view -bT $DB.unzip - | samtools sort -o $Prefix.sorted.bam -
	samtools index $Prefix.sorted.bam
        samtools fastq -f 4 -0 $Prefix.UnMapped.SE.fq $Prefix.sorted.bam
        samtools view $Prefix.sorted.bam | awk -v var="$string" '$3 ~ var' | samtools view -bT $DB.unzip | samtools fastq -0 $Prefix-RFA/$Prefix.Mapped.SE.fq -
	fi
rm $DB.unzip

else
	if [[ -n $R1 && -n $R2 ]]; then
	bwa mem -t $Threads $DB $R1 $R2 | samtools view -bT $DB - | samtools sort -o $Prefix.sorted.bam $Prefix.bam
	samtools index $Prefix.sorted.bam
	samtools fastq -f 4 -1 $Prefix.UnMapped.1.fq -2 $Prefix.UnMapped.2.fq $Prefix.sorted.bam
	samtools view $Prefix.sorted.bam | awk -v var="$string" '$3 ~ var' | samtools view -bT $DB | samtools fastq -1 $Prefix.Mapped.1.fq -2 $Prefix.Mapped.2.fq -
	fi
	if [[ -n $SE ]]; then
	bwa mem -t $Threads $DB $SE | samtools view -bT $DB - | samtools sort -o $Prefix.sorted.bam $Prefix.bam
	samtools index $Prefix.sorted.bam
        samtools fastq -f 4 -0 $Prefix.UnMapped.SE.fq  $Prefix.sorted.bam 
        samtools view $Prefix.sorted.bam | awk -v var="$string" '$3 ~ var' | samtools view -bT $DB | samtools fastq -o $Prefix-RFA/$Prefix.Mapped.SE.fq 
	fi
fi

#repair fastq files
#Need to add subroutine, to only add info when necessary (something like check a read, if it ends \1 then).
pairfq_lite.pl addinfo -i $Prefix.Mapped.1.fq -o $Prefix.Mapped.info.1.fq -p 1
pairfq_lite.pl addinfo -i $Prefix.Mapped.2.fq -o $Prefix.Mapped.info.2.fq -p 2
pairfq_lite.pl addinfo -i $Prefix.UnMapped.1.fq -o $Prefix.UnMapped.info.1.fq -p 1
pairfq_lite.pl addinfo -i $Prefix.UnMapped.2.fq -o $Prefix.UnMapped.info.2.fq -p 2


pairfq_lite.pl makepairs -f $Prefix.Mapped.info.1.fq -r $Prefix.Mapped.info.2.fq -fp $Prefix-RFA/$Prefix.Mapped.RePair.1.fq -rp $Prefix-RFA/$Prefix.Mapped.RePair.2.fq -fs $Prefix-RFA/$Prefix.Mapped.NoPair.1.fq -rs $Prefix-RFA/$Prefix.Mapped.NoPair.2.fq
pairfq_lite.pl makepairs -f $Prefix.UnMapped.info.1.fq -r $Prefix.UnMapped.info.2.fq -fp $Prefix-RFA/$Prefix.UnMapped.RePair.1.fq -rp $Prefix-RFA/$Prefix.UnMapped.RePair.2.fq -fs $Prefix-RFA/$Prefix.UnMapped.NoPair.1.fq -rs $Prefix-RFA/$Prefix.UnMapped.NoPair.2.fq
#Clean up - Saves Gigs of disk space!
#rm $Prefix.UnMapped.NoPair.1.fq
#rm $Prefix.UnMapped.NoPair.2.fq
#rm $Prefix.Mapped.NoPair.1.fq
#rm $Prefix.Mapped.NoPair.2.fq
rm $Prefix.Mapped.1.fq
rm $Prefix.Mapped.2.fq
rm $Prefix.UnMapped.1.fq
rm $Prefix.UnMapped.2.fq
rm $Prefix.sorted.bam
rm $Prefix.sorted.bam.bai
rm $Prefix.Mapped.info.1.fq 
rm $Prefix.Mapped.info.2.fq 
rm $Prefix.UnMapped.info.1.fq 
rm $Prefix.UnMapped.info.2.fq
exit
