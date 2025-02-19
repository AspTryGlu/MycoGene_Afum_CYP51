#!/bin/bash -l


for n in $(cat samples-paths.txt)
do
    sample=$(echo $n | cut -d ":" -f 1)

    
    blastn -query query_seq.fasta -subject denovo_assemblies/$sample/scaffolds.fasta -out blast_$sample -outfmt 6
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

    samtools faidx denovo_assemblies/$sample/scaffolds.fasta $scaff:$start_pos-$end_pos -o seq_$name.fasta
    samtools faidx -i denovo_assemblies/$sample/scaffolds.fasta $scaff:$end_pos-$start_pos -o seq_rev_$name.fasta
done
for file in seq_*;
   do
       sed -i "s/>.*/>${file%%.*}/" "$file"
done
cat seq_* > gene_multifasta.fasta
awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' gene_multifasta.fasta > gene_output.fasta
cat query_seq.fasta gene_output.fasta > gene_output2.fasta
clustalo -i gene_output2.fasta -o aln_gene_output.fasta

rm ../denovo_assemblies/$sample/assembly_graph*
rm ../denovo_assemblies/$sample/K127/assembly_graph*

mkdir intermediate-outputs
mv seq_* intermediate-outputs
mv blast_* intermediate-outputs
mv best_hit_* intermediate-outputs
mv SRR* intermediate-outputs
mv outputfiles-blast intermediate-outputs
mv for_clustal.fasta intermediate-outputs
mv gene_multifasta.fasta intermediate-outputs
mv qc_header intermediate-outputs
