#!/bin/bash -l

mkdir trimmed_reads
mkdir denovo_assemblies

for n in $(cat samples-paths.txt)
do
    sample=$(echo $n | cut -d ":" -f 1)
    path1=$(echo $n | cut -d ":" -f 2)
    path2=$(echo $n | cut -d ":" -f 3)

    FaQCs -1 $path1 -2 $path2  --prefix $sample -d trimmed_reads/
    spades.py -1 trimmed_reads/$sample.1.trimmed.fastq -2 trimmed_reads/$sample.2.trimmed.fastq -k 127  --only-assembler -o denovo_assemblies/$sample
done

for gene in $(cat gene_list)
do
        mkdir $gene
        cd $gene
        cp ../MycoGene_blastn_v1.1.sh .
        cp ../samples-paths.txt .
        cp ../queries/$gene ./query_seq.fasta
        bash MycoGene_blastn_v1.1.sh
        sed -i "s/>seq/>${gene}/" gene_output.fasta
        cd ../

done
