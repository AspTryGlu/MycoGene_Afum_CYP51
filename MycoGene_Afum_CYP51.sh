#!/bin/bash

mkdir trimmed_reads
mkdir denovo_assemblies

for n in $(cat samples-paths.txt)
do
    sample=$(echo $n | cut -d ":" -f 1)
    path1=$(echo $n | cut -d ":" -f 2)
    path2=$(echo $n | cut -d ":" -f 3)

    FaQCs -1 $path1 -2 $path2  --prefix $sample -d trimmed_reads/
    spades.py -1 trimmed_reads/$sample.1.trimmed.fastq -2 trimmed_reads/$sample.2.trimmed.fastq -k 127  --only-assembler -o denovo_assemblies/$sample &
    wait
    blastn -query query_seq.fasta -subject denovo_assemblies/$sample/scaffolds.fasta -out blast_$sample -outfmt 6    
    sed -i 's,^,'"$sample"'\t,' "blast_$sample"
    head -n 1 blast_$sample > best_hit_$sample

    tblastn -query query_aa.fasta -subject denovo_assemblies/$sample/scaffolds.fasta -max_target_seqs 1 -out tblastn_$sample  -outfmt "6 delim=, qstart  sseq"

    sed -i 's/,/\t/g' tblastn_$sample
    sort -k 1 tblastn_$sample |cut -f2 | tr -d '\n' > prot_seq_$sample.fasta
    sed -i 's/EYCFLNRQ//' prot_seq_$sample.fasta
    sed -i 's,^,\>'"$sample"'\n,' "prot_seq_$sample.fasta"

done    

cat best_hit_* > outputfiles-blast
sed -i 's/\t/:/g' outputfiles-blast
echo "sample:qseqid:sseqid:pident:length:mismatch:gapopen:qstart:qend:sstart:send:evalue:bitscore" > qc_header
cat qc_header outputfiles-blast |sed 's/:/\t/g' > QC-report.txt

for line in $(cat outputfiles-blast)
do
    name=$(echo $line | cut -d ":" -f 1)
    scaff=$(echo $line | cut -d ":" -f 3)
    start_pos=$(echo $line | cut -d ":" -f 10)
    end_pos=$(echo $line | cut -d ":" -f 11)
    #upstream length
    new_start_pos=$(($start_pos-500))
    new2_start_pos=$(($start_pos+500))
    
    echo "For sample $name the start POS is $start_pos, the end POS is $end_pos, for including upstream region new start POS is $new_start_pos if gene is reverse New start POS is $new2_start_pos "
    
    samtools faidx ./$name/scaffolds.fasta $scaff:$start_pos-$end_pos -o seq_$name.fasta
    samtools faidx -i ./$name/scaffolds.fasta $scaff:$end_pos-$start_pos -o seq_rev_$name.fasta 

    samtools faidx denovo_assemblies/$name/scaffolds.fasta $scaff:$new_start_pos-$end_pos -o wnoncoding_$name.fasta
    samtools faidx -i denovo_assemblies/$name/scaffolds.fasta $scaff:$end_pos-$new2_start_pos -o wnoncoding_rev_$name.fasta

done
for file in seq_*;
   do
       sed -i "s/>.*/>${file%%.*}/" "$file" 
done
cat seq_* > gene_multifasta.fasta
awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' gene_multifasta.fasta > gene_output.fasta
cat query_seq.fasta gene_output.fasta > gene2_output.fasta
clustalo -i gene2_output.fasta -o aln_gene_output.fasta

for f in prot_seq_*; do cat -- "$f"; printf "\n"; done > protein_output.fasta
cat query_aa.fasta protein_output.fasta > protein2_output.fasta
clustalo -i protein2_output.fasta -o aln_protein_output.fasta


for file in wnoncoding_*;
   do
       sed -i "s/>.*/>${file%%.*}/" "$file"
done
cat wnoncoding_* > wnoncoding_gene_multifasta.fasta
awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' wnoncoding_gene_multifasta.fasta > wnoncoding_gene_output.fasta
cat query_seq.fasta wnoncoding_gene_output.fasta > wnoncoding2_gene_output.fasta
clustalo -i wnoncoding2_gene_output.fasta -o aln_wnoncoding_gene_output.fasta

rm $sample/assembly_graph*
rm $sample/K127/assembly_graph*

mkdir intermediate-outputs
mv seq_* intermediate-outputs
mv blast_* intermediate-outputs
mv best_hit_* intermediate-outputs
mv outputfiles-blast intermediate-outputs
mv gene* intermediate-outputs
mv prot* intermediate-outputs
mv gene_multifasta.fasta intermediate-outputs
mv qc_header intermediate-outputs
mv wnoncoding* intermediate-outputs
mv tblastn_* intermediate-outputs

#vizualization mView
#./mview -in fasta -html head -css on -coloring mismatch -colormap red aln_protein_output.fasta > data_protein_aln.html
./mview -in fasta -html head -css on -coloring mismatch -colormap red -find 'NGKL|DVVY|MMIT|ISYG|IKYG|GFTP|PINFM|EVVDY|LPFG' aln_protein_output.fasta > data2_protein_aln.html
./mview -in fasta -html head -css on -coloring identity aln_gene_output.fasta > data_gene_aln.html
./mview -in fasta -html head -css on -coloring identity aln_wnoncoding_gene_output.fasta > data_upstream-CPY51_aln.html
