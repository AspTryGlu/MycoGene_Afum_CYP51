from Bio import SeqIO

def gap_length_between_patterns(sequence, pattern1, pattern2):
    """Find the length of the gap between two specific patterns in the sequence."""
    index1 = sequence.find(pattern1)
    index2 = sequence.find(pattern2, index1 + len(pattern1))  # Find pattern2 after pattern1

    # Initialize gap length
    gap_length = None

    if index1 != -1 and index2 != -1:
        # Count the gap length between the end of pattern1 and the start of pattern2
        gap_length = 0
        for char in sequence[index1 + len(pattern1):index2]:
            if char == '-':
                gap_length += 1
    else:
        if index1 == -1:
            print(f"Pattern 1 not found in sequence: {sequence[:50]}...")  # Debug: Print part of the sequence
        if index2 == -1:
            print(f"Pattern 2 not found in sequence: {sequence[:50]}...")  # Debug: Print part of the sequence

    return gap_length  # Return gap length (None if patterns not found)

def find_gap_between_patterns(fasta_file, pattern1, pattern2):
    """Find the gap length between two specific patterns in all sequences."""
    results = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        gap_length = gap_length_between_patterns(sequence, pattern1, pattern2)
        results.append((record.id, gap_length))

    return results

def determine_tr(gap_length):
    """Determine the TR value based on gap length."""
    if gap_length is None:
        return "no-TR"  # Return the message directly
    elif gap_length == 0:
        return "TR53"
    elif gap_length == 7:
        return "TR46"
    elif gap_length == 19:
        return "TR34"
    elif gap_length == 53:
        return "no-TR"  # Ensure this condition is checked
    else:
        return "Unknown"  # For any other gap lengths

def print_results(results):
    """Print the results."""
    print(f"{'ID':<20} {'Gap Length':<30} {'TR'}")
    print("="*70)
    for entry in results:
        seq_id, gap_length = entry
        tr_value = determine_tr(gap_length)
        print(f"{seq_id:<20} {gap_length if gap_length is not None else 'N/A':<30} {tr_value}")

def save_results_to_file(results, output_file):
    """Save the results to a specified output file."""
    with open(output_file, 'w') as f:
        f.write(f"{'ID':<20} {'Gap Length':<30} {'TR'}\n")
        f.write("="*70 + "\n")
        for entry in results:
            seq_id, gap_length = entry
            tr_value = determine_tr(gap_length)
            f.write(f"{seq_id:<20} {gap_length if gap_length is not None else 'N/A':<30} {tr_value}\n")

# Usage
fasta_file = "aln_wnoncoding_gene_output.fasta"  # Replace with your FASTA file path
output_file = "Tandem_repeats_results_from_aln.txt"
pattern1 = "ATAATCGCAGCACCACTTCAGAGTTGTCTAGAATCACGCGGTCCGGATGT"  # First pattern
pattern2 = "AGTTGCCTAATTACTAAGGTGTAGTTCCAGCATACCATACACCCT"  # Second pattern
results = find_gap_between_patterns(fasta_file, pattern1, pattern2)

print_results(results)
save_results_to_file(results, output_file)
