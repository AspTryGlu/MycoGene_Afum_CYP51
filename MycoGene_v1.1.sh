#!/bin/bash


for n in $(cat samples-paths.txt)
do
    sample=$(echo $n | cut -d ":" -f 1)
    path1=$(echo $n | cut -d ":" -f 2)
    path2=$(echo $n | cut -d ":" -f 3)

    spades.py -1 $path1 -2 $path2 -k 127  --only-assembler -o $sample
    tblastn -query query_seq.fasta -subject $sample/scaffolds.fasta -out blast_$sample -outfmt 6    
    sed -i 's,^,'"$sample"'\t,' "blast_$sample"
    head -n 1 blast_$sample > best_hit_$sample
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

    samtools faidx ./$name/scaffolds.fasta $scaff:$start_pos-$end_pos -o seq_$name.fasta
    samtools faidx -i ./$name/scaffolds.fasta $scaff:$end_pos-$start_pos -o seq_rev_$name.fasta    
done
for file in seq_*;
   do
       sed -i "s/>.*/>${file%%.*}/" "$file" 
done
cat seq_* > gene_multifasta.fasta
awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' gene_multifasta.fasta > gene_output.fasta
clustalo -i gene_output.fasta -o aln_gene_output.fasta

transeq -sequence gene_output.fasta -outseq protein_output.fasta
cat query_seq.fasta protein_output.fasta > for_clustal.fasta
clustalo -i for_clustal.fasta -o aln_protein_output.fasta

rm $sample/assembly_graph*
rm $sample/K127/assembly_graph*

mkdir intermediate-outputs
mv seq_* intermediate-outputs
mv blast_* intermediate-outputs
mv best_hit_* intermediate-outputs
mv SRR* intermediate-outputs
mv outputfiles-blast intermediate-outputs
mv for_clustal.fasta intermediate-outputs
mv gene_multifasta.fasta intermediate-outputs
mv qc_header intermediate-outputs
