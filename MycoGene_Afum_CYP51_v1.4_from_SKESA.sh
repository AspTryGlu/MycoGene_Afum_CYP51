#!/bin/bash -l

mkdir trimmed_reads
mkdir denovo_assemblies

SPADES_IMAGE_FILE="spades_3.15.5.sif"
SPADES_IMAGE_URL="docker://quay.io/staphb/spades:3.15.5"
if [ -f "$SPADES_IMAGE_FILE" ]; then
    echo "Image $SPADES_IMAGE_FILE already exists. Skipping pull."
else
    echo "Image $SPADES_IMAGE_FILE does not exist. Pulling the image..."
    singularity pull "$SPADES_IMAGE_URL"
fi
SPADES_THREADS=$(( $(nproc) - 2 ))

CLUSTAL_IMAGE_FILE="clustalo_1.2.4.sif"
CLUSTAL_IMAGE_URL="docker://quay.io/staphb/clustalo:1.2.4"
if [ -f "$CLUSTAL_IMAGE_FILE" ]; then
    echo "Image $CLUSTAL_IMAGE_FILE already exists. Skipping pull."
else
    echo "Image $CLUSTAL_IMAGE_FILE does not exist. Pulling the image..."
    singularity pull "$CLUSTAL_IMAGE_URL"
fi


for n in $(cat samples-paths.txt)
do
    sample=$(echo $n | cut -d ":" -f 1)
    path1=$(echo $n | cut -d ":" -f 2)
    path2=$(echo $n | cut -d ":" -f 3)

    #FaQCs -1 $path1 -2 $path2  --prefix $sample -d trimmed_reads/
    #singularity exec "$SPADES_IMAGE_FILE" spades.py -1 trimmed_reads/$sample.1.trimmed.fastq -2 trimmed_reads/$sample.2.trimmed.fastq -k 127  --only-assembler -o denovo_assemblies/$sample -t "$SPADES_THREADS" &
    #wait
    blastn -query query_seq.fasta -subject denovo_assemblies/$sample.fa -out blast_$sample -outfmt 6
    sed -i 's,^,'"$sample"'\t,' "blast_$sample"
    head -n 1 blast_$sample > best_hit_$sample

    tblastn -query query_aa.fasta -subject denovo_assemblies/$sample.fa -max_target_seqs 1 -out tblastn_$sample  -outfmt "6 delim=, qstart  sseq"

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
    
    samtools faidx denovo_assemblies/$name.fa $scaff:$start_pos-$end_pos -o seq_$name.fasta
    samtools faidx -i denovo_assemblies/$name.fa $scaff:$end_pos-$start_pos -o seq_rev_$name.fasta 

    samtools faidx denovo_assemblies/$name.fa $scaff:$new_start_pos-$end_pos -o wnoncoding_$name.fasta
    samtools faidx -i denovo_assemblies/$name.fa $scaff:$end_pos-$new2_start_pos -o wnoncoding_rev_$name.fasta

done
for file in seq_*;
   do
       sed -i "s/>.*/>${file%%.*}/" "$file" 
done
cat seq_* > gene_multifasta.fasta
awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' gene_multifasta.fasta > gene_output.fasta
cat query_seq.fasta gene_output.fasta > gene2_output.fasta
#clustalo -i gene2_output.fasta -o aln_gene_output.fasta
#clustalw2 -INFILE=gene2_output.fasta -align -OUTFILE=aln_gene_output.fasta -OUTPUT=FASTA
singularity exec "$CLUSTAL_IMAGE_FILE" clustalo -i gene2_output.fasta --dealign -o aln_gene_output.fasta --force

for f in prot_seq_*; do cat -- "$f"; printf "\n"; done > protein_output.fasta
cat query_aa.fasta protein_output.fasta > protein2_output.fasta
#clustalo -i protein2_output.fasta -o aln_protein_output.fasta
#clustalw2 -INFILE=protein2_output.fasta -align -OUTFILE=aln_protein_output.fasta -OUTPUT=FASTA
singularity exec "$CLUSTAL_IMAGE_FILE" clustalo -i protein2_output.fasta --dealign -o aln_protein_output.fasta --force

for file in wnoncoding_*;
   do
       sed -i "s/>.*/>${file%%.*}/" "$file"
done
cat wnoncoding_* > wnoncoding_gene_multifasta.fasta
awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' wnoncoding_gene_multifasta.fasta > wnoncoding_gene_output.fasta
cat ctrl_TR wnoncoding_gene_output.fasta > wnoncoding2_gene_output.fasta
#clustalo -i wnoncoding2_gene_output.fasta -o aln_wnoncoding_gene_output.fasta
#clustalw2 -INFILE=wnoncoding2_gene_output.fasta -align -OUTFILE=aln_wnoncoding_gene_output.fasta -OUTPUT=FASTA
singularity exec "$CLUSTAL_IMAGE_FILE" clustalo -i wnoncoding2_gene_output.fasta --dealign  -o aln_wnoncoding_gene_output.fasta --force

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

##table output


module load python/3.11.3
pip install biopython


python TR_table_v7.py
python Mismatches_table_v3.py

##Table TR

blastn -query query_long_TR53 -subject intermediate-outputs/wnoncoding2_gene_output.fasta -gapopen 2 -gapextend 1 -out blastn_gapopenmod_TR_determination -outfmt 7

awk '{
    col2 = $2;  # Store column 2
    col9 = $9;  # Store column 9
    # Apply transformations based on the value of column 9
    if (col9 == 1) {
        col9 = "TR53";
    } else if (col9 == 11) {
        col9 = "TR46";
    } else if (col9 == 14) {
        col9 = "TR34";
    } else if (col9 == 78) {
        col9 = "nonTR";
    }
    # Print column 2 and the transformed column 9
    print col2, col9;
}' blastn_gapopenmod_TR_determination > TR_results.txt
