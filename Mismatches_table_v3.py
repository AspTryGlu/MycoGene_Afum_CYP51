#DISCLAIMER: This script includes code that has been partially or fully generated using AI tools. While efforts have been made to ensure accuracy and functionality, users are advised to review and test the code thoroughly before implementing it.

from Bio import SeqIO

def find_mismatches(fasta_file, reference_id):
    """Find mismatches between sequences and a reference sequence, reporting ungapped positions."""
    mismatches = []
    reference_sequence = None

    # First, retrieve the reference sequence
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id == reference_id:
            reference_sequence = str(record.seq)
            break

    if reference_sequence is None:
        raise ValueError(f"Reference sequence with ID '{reference_id}' not found.")

    # Prepare an ungapped reference sequence and its positions
    ungapped_reference = []
    ungapped_positions = []
    for pos, aa in enumerate(reference_sequence, start=1):
        if aa != '-':  # Ignore gaps in the reference
            ungapped_reference.append(aa)
            ungapped_positions.append(pos)

    # Now, compare each sequence to the ungapped reference sequence
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id != reference_id:  # Skip the reference itself
            sequence = str(record.seq)
            ungapped_index = 0  # Index for ungapped reference positions
            sequence_index = 0  # Index for the sequence

            for ref_aa in ungapped_reference:
                # Move the sequence index to the next non-gap amino acid
                while sequence_index < len(sequence) and sequence[sequence_index] == '-':
                    sequence_index += 1

                # Check if we are still within the bounds of the sequence
                if sequence_index < len(sequence):
                    seq_aa = sequence[sequence_index]
                    if ref_aa != seq_aa:
                        # Report mismatch with position from ungapped reference
                        mismatch_position = ungapped_positions[len(ungapped_reference) - len(ungapped_reference) + ungapped_index]
                        mismatches.append((record.id, f"{ref_aa}{mismatch_position}{seq_aa}"))
                
                # Move to the next position in the ungapped reference
                ungapped_index += 1
                sequence_index += 1  # Always move to the next position in the sequence

    return mismatches

def save_mismatch_results_to_file(mismatches, output_file):
    """Save the mismatch results to a specified output file."""
    with open(output_file, 'w') as f:
        f.write(f"{'Sequence ID':<20} {'Mismatch'}\n")
        f.write("="*50 + "\n")
        for entry in mismatches:
            seq_id, mismatch = entry
            f.write(f"{seq_id:<20} {mismatch}\n")

# Example usage
fasta_file = "aln_protein_output.fasta"  # Replace with your protein FASTA file path
reference_id = "AAK73659"  # Reference sequence ID
mismatch_results_file = "mismatch_results.txt"  # Specify the output file name for mismatches

# Find mismatches and save to file
mismatches = find_mismatches(fasta_file, reference_id)
save_mismatch_results_to_file(mismatches, mismatch_results_file)

print(f"Mismatches saved to {mismatch_results_file}")
