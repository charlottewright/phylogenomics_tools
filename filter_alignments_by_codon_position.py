#!/usr/bin/env python3
#%%
import argparse

def read_fasta_file(file_path):
    sequences = {}  # Dictionary to store the pairs of sequence IDs and sequences
    current_id = None  # To store the current ID
    current_sequence = []  # To accumulate the sequence

    valid_bases = {'A', 'T', 'C', 'G', 'a', 't', 'c', 'g', 'r', 'y', 'R','Y'}  # Set of valid bases

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()  # Remove any leading/trailing whitespace
            if line.startswith('>'):  # If line starts with '>', it's an ID
                if current_id is not None: # If we already have an ID, save the previous sequence
                    sequences[current_id] = ''.join(current_sequence)
                # Set the new ID, and reset the current sequence
                current_id = line[1:]  # Remove the '>' character
                current_sequence = []
            else:
                # Keep only valid bases (A, T, C, G) and ignore '-'
                filtered_sequence = [char for char in line if char in valid_bases]
                current_sequence.extend(filtered_sequence)

        # Save the last sequence after the loop finishes
        if current_id is not None:
            sequences[current_id] = ''.join(current_sequence)
    return sequences

def keep_every_n_letters(sequences, n): # keep every third letter (i.e a codon position), given a starting number e.g. 1 means keep every first codon position
    modified_sequences = {}  # Dictionary to store ID and the modified sequence

    for seq_id, sequence in sequences.items():
        start = n - 1  # Adjust n to be zero-indexed (1 -> 0, 2 -> 1, etc.)
        modified_sequence = sequence[start::3]  # Keep every third letter, starting from the n-th position
        modified_sequences[seq_id] = modified_sequence
    return modified_sequences

def interleave_sequences(dict1, dict2):
    merged_sequences = {}

    # Iterate over the sequence IDs (keys) in dict1
    for seq_id in dict1:
        seq1 = dict1[seq_id]
        seq2 = dict2[seq_id]

        merged_sequence = ''.join(a + b for a, b in zip(seq1, seq2)) # Interleave the letters from both sequences
        merged_sequences[seq_id] = merged_sequence # Store the new sequence in the merged_sequences dictionary

    return merged_sequences

if __name__ == "__main__":
    SCRIPT = "test_for_base_composition.py"
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', 
                        type=str, 
                        required=True, 
                        help="Path to the input file")

    # Add the output prefix argument
    parser.add_argument('-p', '--prefix', 
                        type=str, 
                        required=True, 
                        help="Prefix for the output files")

    parser.add_argument('-P', '--positions', 
                        type=str, 
                        required=True, 
                        help="Prefix for the output files")
    # Parse the command-line arguments
    args = parser.parse_args()

    input_file = args.file
    output_prefix = args.prefix
    positions_to_keep = args.positions # parse as list
    positions_to_keep = [int(x) for x in positions_to_keep.split(',')]
    print(f"[+]   Reading input alignment file: {input_file}")
    print(f"[+]   Filtering alignments to only keep codon positions: {positions_to_keep}")

    sequences_dict = read_fasta_file(input_file)
    sequences_dict = {seq_id: sequence.lower() for seq_id, sequence in sequences_dict.items()}

    if len(positions_to_keep) == 1: # i.e if one one position is required, this is simple - just extract that position
        sequences_filtered = keep_every_n_letters(sequences_dict, positions_to_keep[0]) # keep just the codon position which isnt biased
    elif set(positions_to_keep) == set([1,2]):
        sequences_pos1 = keep_every_n_letters(sequences_dict, 1) # keep codon pos 1
        sequences_pos2 = keep_every_n_letters(sequences_dict, 2) # keep codon pos 2
        sequences_filtered = interleave_sequences(sequences_pos1, sequences_pos2) # order matters, start with first codon positions then second each time
    elif set(positions_to_keep) == set([1,3]):
        sequences_pos1 = keep_every_n_letters(sequences_dict, 1) # keep codon pos 1
        sequences_pos2 = keep_every_n_letters(sequences_dict, 3) # keep codon pos 3
        sequences_filtered = interleave_sequences(sequences_pos1, sequences_pos2)  
    elif set(positions_to_keep) == set([2,3]):
        sequences_pos1 = keep_every_n_letters(sequences_dict, 2) # keep codon pos 2
        sequences_pos2 = keep_every_n_letters(sequences_dict, 3) # keep codon pos 3
        sequences_filtered = interleave_sequences(sequences_pos1, sequences_pos2) 

    output_file = output_prefix + '.fa'
    print(f"[+]   Writing filtered alignments to: {output_file}")
    with open(output_file, 'w') as output:
        for seq_id, sequence in sequences_filtered.items():
            # Write each sequence in the desired format
            output.write(f">{seq_id}\n{sequence}\n")        
