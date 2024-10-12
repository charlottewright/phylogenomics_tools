#!/usr/bin/env python3
#%%
import scipy.stats as stats
import argparse
import sys

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

def count_bases(sequences):
    base_counts = {}  # Dictionary to store ID and the base counts as a tuple
    for seq_id, sequence in sequences.items(): # Count the occurrences of each base
        A_count = sequence.count('a')
        C_count = sequence.count('c')
        T_count = sequence.count('t')
        G_count = sequence.count('g')
        R_count = sequence.count('r')
        Y_count = sequence.count('y')
        # Store the counts as a tuple in the dictionary
        base_counts[seq_id] = (A_count, C_count, T_count, G_count, R_count, Y_count)
    return base_counts


def calculate_proportions(sequences):
    base_proportions = {}  # Dictionary to store ID and the base proportions as a tuple

    for seq_id, sequence in sequences.items():
        total_bases = len(sequence)  # Total length of the sequence (only valid bases)
        
        if total_bases == 0:
            # If no valid bases, set proportions to zero
            base_proportions[seq_id] = (0, 0, 0, 0)
            continue

        # Count the occurrences of each base
        A_count = sequence.count('a')
        C_count = sequence.count('c')
        T_count = sequence.count('t')
        G_count = sequence.count('g')
        R_count = sequence.count('r')
        Y_count = sequence.count('y')

        # Calculate proportions for each base
        A_prop = A_count / total_bases
        C_prop = C_count / total_bases
        T_prop = T_count / total_bases
        G_prop = G_count / total_bases
        R_prop = R_count / total_bases
        Y_prop = Y_count / total_bases
        # Store the proportions as a tuple in the dictionary
        base_proportions[seq_id] = (A_prop, C_prop, T_prop, G_prop, R_prop, Y_prop)

    return base_proportions


def calculate_average_proportions_by_counts(sequences):
    total_A, total_C, total_T, total_G, total_R, total_Y = 0, 0, 0, 0, 0, 0
    total_bases = 0  # Total number of valid bases (A, C, T, G)

    for sequence in sequences.values():
        # Count the occurrences of each base
        A_count = sequence.count('a')
        C_count = sequence.count('c')
        T_count = sequence.count('t')
        G_count = sequence.count('g')
        R_count = sequence.count('r')
        Y_count = sequence.count('y')
        # Sum the counts across all sequences
        total_A += A_count
        total_C += C_count
        total_T += T_count
        total_G += G_count
        total_R += R_count
        total_Y += Y_count
        # Update the total number of valid bases
        total_bases += (A_count + C_count + T_count + G_count + R_count + Y_count)

        # Calculate the proportion of each base
        A_prop = total_A / total_bases
        C_prop = total_C / total_bases
        T_prop = total_T / total_bases
        G_prop = total_G / total_bases
        R_prop = total_R / total_bases
        Y_prop = total_Y / total_bases
    return A_prop, C_prop, T_prop, G_prop, R_prop, Y_prop

def keep_every_n_letters(sequences, n): # keep every third letter (i.e a codon position), given a starting number e.g. 1 means keep every first codon position
    modified_sequences = {}  # Dictionary to store ID and the modified sequence

    for seq_id, sequence in sequences.items():
        start = n - 1  # Adjust n to be zero-indexed (1 -> 0, 2 -> 1, etc.)
        modified_sequence = sequence[start::3]  # Keep every third letter, starting from the n-th position
        modified_sequences[seq_id] = modified_sequence
    return modified_sequences

def calculate_expected_counts(sequences, avg_proportions): # calculate the expected number of A, C, T, G in each sequence based on the average proportions of A,C,T,G across all sequences
    expected_counts = {}  # Dictionary to store expected counts for each sequence

    avg_A, avg_C, avg_T, avg_G, avg_R, avg_Y = avg_proportions  # Unpack average proportions

    for seq_id, sequence in sequences.items():
        sequence_length = len(sequence)  # Length of the sequence

        # Calculate expected counts for A, C, T, G
        expected_A = sequence_length * avg_A
        expected_C = sequence_length * avg_C
        expected_T = sequence_length * avg_T
        expected_G = sequence_length * avg_G
        expected_R = sequence_length * avg_R
        expected_Y = sequence_length * avg_Y
        # Store the expected counts as a tuple in the dictionary
        expected_counts[seq_id] = (expected_A, expected_C, expected_T, expected_G, expected_R, expected_Y)
    return expected_counts

def calculate_chi2_across_all_sequences_of_gene(sequence_dict, observed_dict, expected_dict): 
    # Initialize lists to store all observed and expected values across all sequences and bases
    all_observed, all_expected = [], []

    # Flatten observed and expected counts into single lists
    all_letter_types = ''.join(value for value in sequence_dict.values()) # lets find all letters across all sequences
    unique_letters = set(all_letter_types) # get unique letters
    letter_to_index_dict = {'a':0, 't':1, 'c':2, 'g':3, 'r':4, 'y':5}
    for seq_id in observed_dict:
        for letter in unique_letters: # iterate over A,C,T,G etc
            index_pos = letter_to_index_dict[letter]
            all_observed.append(observed_dict[seq_id][index_pos])
            all_expected.append(expected_dict[seq_id][index_pos])
    # Perform the chi-squared test on the flattened data
    chi2_stat, p_value = stats.chisquare(f_obs=all_observed, f_exp=all_expected)
    all_letter_types = ''.join(value for value in sequence_dict.values()) # lets find all letters across all sequences
    unique_letters = set(all_letter_types) # get unique letters
    degrees_of_freedom = (len(observed_dict) -1 )*(len(unique_letters)-1)

    return(chi2_stat, p_value, degrees_of_freedom)


def read_fasta_file(file_path):
    sequences = {}  # Dictionary to store ID and sequence pairs
    current_id = None  # To store the current ID
    current_sequence = []  # To accumulate the sequence

    valid_bases = {'A', 'T', 'C', 'G', 'a', 't', 'c', 'g'}  # Set of valid bases

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()  # Remove any leading/trailing whitespace
            if line.startswith('>'):  # If line starts with '>', it's an ID
                if current_id is not None:
                    # If we already have an ID, save the previous sequence
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


# Function to replace every n-th character in the sequence
def replace_nth_letter(sequence, n):
    modified_sequence = list(sequence)  # Convert to a list for easier manipulation
    for i in range(n, len(modified_sequence), 3):  # Start from n and step by 3
        if modified_sequence[i] == 'a':
            modified_sequence[i] = 'r'
        elif modified_sequence[i] == 'g':
            modified_sequence[i] = 'r'
        elif modified_sequence[i] == 't':
            modified_sequence[i] = 'y'
        elif modified_sequence[i] == 'c':
            modified_sequence[i] = 'y'
    return ''.join(modified_sequence)  # Join back into a string

def interleave_sequences(dict1, dict2):
    merged_sequences = {}

    # Iterate over the sequence IDs (keys) in dict1
    for seq_id in dict1:
        seq1 = dict1[seq_id]
        seq2 = dict2[seq_id]

        merged_sequence = ''.join(a + b for a, b in zip(seq1, seq2)) # Interleave the letters from both sequences
        merged_sequences[seq_id] = merged_sequence # Store the new sequence in the merged_sequences dictionary

    return merged_sequences


def write_partition_file(sequences_dict, num_pos, file_name):
    # Open the file in write mode
    with open(file_name, 'w') as file:
        # Iterate through the dictionary of sequence IDs and sequences
        for seq_id, sequence in sequences_dict.items():
            # Calculate the length of the current sequence
            length = len(sequence)
            # Write lines for this sequence ID
            for i in range(1, num_pos + 1):
                line = f"DNA, {seq_id}_part{i} = {i}-{length}/{i},\n"
                file.write(line)
    return()

#%%
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
    parser.add_argument('-r','--recode',
                        help='Enable recoding of codon positions which show composition bias to RY',
                        type=str,
                        required=False,
                        default=False)
    parser.add_argument('-d','--delete',
                        help='Enable removal of codon positions which show composition bias after RY recoding',
                        type=str,
                        required=False,
                        default=False)
    parser.add_argument('-P','--partition',
                        help='Generate partition files for the input alignment and, if generated, for the alignment which has one or more codon positions removed.',
                        type=str,
                        required=False,
                        default=False)
    # Parse the command-line arguments
    args = parser.parse_args()

    input_file = args.file
    output_prefix = args.prefix
    recode_requested = args.recode
    remove_requested = args.delete
    partitions_requested = args.partition
    # input_file = '../Analysis/phylogeny_redo/backtrans_nucleotide_alignments_trimmed/999at7088.faa.reformatted.aln.trimal'  # Replace with your file path

    print(f"[+]   Reading input alignment file: {input_file}")
    print(f"[+]   Output prefix: {output_prefix}")
    print(f"[+]   RY-recode requested if needed: {remove_requested}")
    print(f"[+]   Removal of codon positions if bias remains after recoding: {remove_requested}")
    print(f"[+]   Partition file(s) requested: {remove_requested}")

    print('')
    sequences_dict = read_fasta_file(input_file)
    sequences_dict = {seq_id: sequence.lower() for seq_id, sequence in sequences_dict.items()}

    # To print the dictionary:
    #for seq_id, sequence in sequences_dict.items():
    #    print(f"ID: {seq_id}, Sequence: {sequence}")

    if partitions_requested != False:
        parititon_file = output_prefix + '.partitions' 
        print(f"[+]   Writing partition file for the alignment to {parititon_file}")
        write_partition_file(sequences_dict, 3, parititon_file)
    # %%
    proportions_dict = calculate_proportions(sequences_dict)
    average_proportions = calculate_average_proportions_by_counts(sequences_dict) # calculate the average proportions of each base across all sequences

    first_codon_sequences_dict = keep_every_n_letters(sequences_dict, 1) # keep just first codon position
    second_codon_sequences_dict = keep_every_n_letters(sequences_dict, 2) # keep just second codon position
    third_codon_sequences_dict = keep_every_n_letters(sequences_dict, 3) # keep just third codon position

    proportions_dict = calculate_proportions(sequences_dict) # across all codon positions
    average_proportions = calculate_average_proportions_by_counts(sequences_dict)
    first_codon_proportions_dict = calculate_proportions(first_codon_sequences_dict)
    first_codon_average_proportions = calculate_average_proportions_by_counts(first_codon_sequences_dict)
    second_codon_proportions_dict = calculate_proportions(second_codon_sequences_dict)
    second_codon_average_proportions = calculate_average_proportions_by_counts(second_codon_sequences_dict)
    third_codon_proportions_dict = calculate_proportions(third_codon_sequences_dict)
    third_codon_average_proportions = calculate_average_proportions_by_counts(third_codon_sequences_dict)

    # now lets calculate expected proportions 
    expected_counts_dict = calculate_expected_counts(sequences_dict, average_proportions)
    first_codon_expected_counts_dict = calculate_expected_counts(first_codon_sequences_dict, first_codon_average_proportions)
    second_codon_expected_counts_dict = calculate_expected_counts(second_codon_sequences_dict, second_codon_average_proportions)
    third_codon_expected_counts_dict = calculate_expected_counts(third_codon_sequences_dict, third_codon_average_proportions)

    base_counts_dict = count_bases(sequences_dict)
    first_codon_counts_dict = count_bases(first_codon_sequences_dict)
    second_codon_counts_dict = count_bases(second_codon_sequences_dict)
    third_codon_counts_dict = count_bases(third_codon_sequences_dict)

    # To print the base counts dictionary:
    #for seq_id, counts in third_codon_counts_dict.items():
    #    print(f"ID: {seq_id}, A: {counts[0]}, C: {counts[1]}, T: {counts[2]}, G: {counts[3]}")
    # %%
    all_bases_chi2_stat, all_bases_p_value, all_bases_df = calculate_chi2_across_all_sequences_of_gene(sequences_dict, base_counts_dict,expected_counts_dict)
    first_codon_chi2_stat, first_codon_p_value, first_codon_df = calculate_chi2_across_all_sequences_of_gene(first_codon_sequences_dict, first_codon_counts_dict,first_codon_expected_counts_dict)
    second_codon_chi2_stat, second_codon_p_value, second_codon_df = calculate_chi2_across_all_sequences_of_gene(second_codon_sequences_dict, second_codon_counts_dict,second_codon_expected_counts_dict)
    third_codon_chi2_stat, third_codon_p_value, third_codon_df = calculate_chi2_across_all_sequences_of_gene(third_codon_sequences_dict, third_codon_counts_dict,third_codon_expected_counts_dict)
    print("[+]   Testing for base composition bias:")
    print('')
    print('Across all codon position:')
    print(f"Chi-squared statistic: {all_bases_chi2_stat} df=({all_bases_df}) P-value: {all_bases_p_value}")

    print('First codon position:')
    print(f"Chi-squared statistic: {first_codon_chi2_stat} df=({first_codon_df}) P-value: {first_codon_p_value}")

    print('Second codon position:')
    print(f"Chi-squared statistic: {second_codon_chi2_stat} df=({second_codon_df}) P-value: {second_codon_p_value}")

    print('Third codon position:')
    print(f"Chi-squared statistic: {third_codon_chi2_stat} df=({third_codon_df}) P-value: {third_codon_p_value}")

    output_file = output_prefix + '.base_composition_results.tsv'
    with open(output_file, 'w') as output:
        output.write('positions\tchi2_stat\tchr_sq_df\tp_value\n')  # header    
        output.write(f'input_alignments\tall_positions\t{all_bases_chi2_stat}\t{all_bases_df}\t{all_bases_p_value}\n') 
        output.write(f'input_alignments\tfirst_codon\t{first_codon_chi2_stat}\t{first_codon_df}\t{first_codon_p_value}\n') 
        output.write(f'input_alignments\tsecond_codon\t{second_codon_chi2_stat}\t{second_codon_df}\t{second_codon_p_value}\n') 
        output.write(f'input_alignments\tthird_codon\t{third_codon_chi2_stat}\t{third_codon_df}\t{third_codon_p_value}\n') 
        print('')
        print(f"[+]   Base composition results written to {output_file}.")

# %%
    output_alignments_file = output_prefix + '.recoded'
    codon_positions_actual = []
    if first_codon_p_value < 0.05: # if codon bias at first codon position
        codon_positions_actual.append(1)
    if second_codon_p_value < 0.05:
        codon_positions_actual.append(2)
    if third_codon_p_value < 0.05:
        codon_positions_actual.append(3)
    codon_positions = [num - 1 for num in codon_positions_actual] # need to minus one from each to account for how numbers are used in python
    # input_file = '../Analysis/phylogeny_redo/backtrans_nucleotide_alignments_trimmed/999at7088.faa.reformatted.aln.trimal'  # Replace with your file path
    if recode_requested == False:
         print('[+] Either RY-recoding was not enabled or no codon positions show composition bias. Exiting script.')
         sys.exit()
    elif len(codon_positions)  == 0: # if no codon positions show bias, we don't need to recode
        print('[+] No codon positions show composition bias. Exiting script.')
        with open(output_file, 'a') as output:
            output.write('Outcome: No_codon_composition_bias.\n') 
        sys.exit()
    else:
        print(f"[+]   Performing RY recoding of codon(s) {codon_positions_actual} due to base composition bias (p < 0.05)")
        # Iterate over the dictionary and replace characters
        modified_sequences_dict = {}
        # Set n to the desired interval
        n = 3  # Replace every 3rd letter
        # Iterate over the original dictionary and replace characters
        recoded_sequences_dict = sequences_dict
        for seq_id, sequence in sequences_dict.items():
            for position in codon_positions:
                modified_sequence = replace_nth_letter(sequence, position)  # Replace every position (first, second etc)
                sequence = modified_sequence
            recoded_sequences_dict[seq_id] = sequence
                #   modified_sequences_dict[seq_id] = modified_sequence.lower()  # Convert to lowercase

        with open(output_alignments_file, 'w') as output:
            for seq_id, sequence in recoded_sequences_dict.items():
                # Write each sequence in the desired format
                output.write(f">{seq_id}\n{sequence}\n")

        print(f"[+]   RY recoded sequences written to {output_alignments_file}.")

#%%
        print('[+]   Repeating calculations of base composition bias on RY recoded sequences.')
        print('')
        proportions_dict = calculate_proportions(recoded_sequences_dict)
        average_proportions = calculate_average_proportions_by_counts(recoded_sequences_dict) # calculate the average proportions of each base across all sequences

        first_codon_sequences_dict = keep_every_n_letters(recoded_sequences_dict, 1) # keep just first codon position
        second_codon_sequences_dict = keep_every_n_letters(recoded_sequences_dict, 2) # keep just second codon position
        third_codon_sequences_dict = keep_every_n_letters(recoded_sequences_dict, 3) # keep just third codon position

        proportions_dict = calculate_proportions(recoded_sequences_dict) # across all codon positions
        average_proportions = calculate_average_proportions_by_counts(recoded_sequences_dict)
        first_codon_proportions_dict = calculate_proportions(first_codon_sequences_dict)
        first_codon_average_proportions = calculate_average_proportions_by_counts(first_codon_sequences_dict)
        second_codon_proportions_dict = calculate_proportions(second_codon_sequences_dict)
        second_codon_average_proportions = calculate_average_proportions_by_counts(second_codon_sequences_dict)
        third_codon_proportions_dict = calculate_proportions(third_codon_sequences_dict)
        third_codon_average_proportions = calculate_average_proportions_by_counts(third_codon_sequences_dict)
        # now lets calculate expected proportions 
        expected_counts_dict = calculate_expected_counts(recoded_sequences_dict, average_proportions)
        first_codon_expected_counts_dict = calculate_expected_counts(first_codon_sequences_dict, first_codon_average_proportions)
        second_codon_expected_counts_dict = calculate_expected_counts(second_codon_sequences_dict, second_codon_average_proportions)
        third_codon_expected_counts_dict = calculate_expected_counts(third_codon_sequences_dict, third_codon_average_proportions)

        base_counts_dict = count_bases(recoded_sequences_dict)
        first_codon_counts_dict = count_bases(first_codon_sequences_dict)
        second_codon_counts_dict = count_bases(second_codon_sequences_dict)
        third_codon_counts_dict = count_bases(third_codon_sequences_dict)

        all_bases_chi2_stat, all_bases_p_value, all_bases_df = calculate_chi2_across_all_sequences_of_gene(sequences_dict, base_counts_dict,expected_counts_dict)
        first_codon_chi2_stat, first_codon_p_value, first_codon_df = calculate_chi2_across_all_sequences_of_gene(first_codon_sequences_dict, first_codon_counts_dict,first_codon_expected_counts_dict)
        second_codon_chi2_stat, second_codon_p_value, second_codon_df = calculate_chi2_across_all_sequences_of_gene(second_codon_sequences_dict, second_codon_counts_dict,second_codon_expected_counts_dict)
        third_codon_chi2_stat, third_codon_p_value, third_codon_df = calculate_chi2_across_all_sequences_of_gene(third_codon_sequences_dict, third_codon_counts_dict,third_codon_expected_counts_dict)

        print('Across all codon position:')
        print(f"Chi-squared statistic: {all_bases_chi2_stat} df=({all_bases_df}) P-value: {all_bases_p_value}")

        print('First codon position:')
        print(f"Chi-squared statistic: {first_codon_chi2_stat} df=({first_codon_df}) P-value: {first_codon_p_value}")

        print('Second codon position:')
        print(f"Chi-squared statistic: {second_codon_chi2_stat} df=({second_codon_df}) P-value: {second_codon_p_value}")

        print('Third codon position:')
        print(f"Chi-squared statistic: {third_codon_chi2_stat} df=({third_codon_df}) P-value: {third_codon_p_value}")

       # output_file = output_prefix + '.base_composition.RY_recoded.tsv' 
        output_file = output_prefix + '.base_composition_results.tsv'
        with open(output_file, 'a') as output:
 #           output.write('positions\tchi2_stat\tchr_sq_df\tp_value\n')  # header    
            output.write(f'RY_recoded_alignments\tall_positions\t{all_bases_chi2_stat}\t{all_bases_df}\t{all_bases_p_value}\n') 
            output.write(f'RY_recoded_alignments\tfirst_codon\t{first_codon_chi2_stat}\t{first_codon_df}\t{first_codon_p_value}\n') 
            output.write(f'RY_recoded_alignments\tsecond_codon\t{second_codon_chi2_stat}\t{second_codon_df}\t{second_codon_p_value}\n') 
            output.write(f'RY_recoded_alignments\tthird_codon\t{third_codon_chi2_stat}\t{third_codon_df}\t{third_codon_p_value}\n') 
            if (first_codon_p_value < 0.05) & (second_codon_p_value < 0.05) & (third_codon_p_value < 0.05):
                output.write('Warning! All three codon positions still show composition bias, even after RY recoding.\n')
            print('')
            print(f"[+]   Base composition results on RY-recoded sequences written to written to {output_file}.")

#%%
        # if after doing RY-recoding, we still have biased positions, we need to remove them.
        codon_positions_actual = []
        if first_codon_p_value < 0.05: # if codon bias at first codon position
            codon_positions_actual.append(1)
        if second_codon_p_value < 0.05:
            codon_positions_actual.append(2)
        if third_codon_p_value < 0.05:
            codon_positions_actual.append(3)
        # input_file = '../Analysis/phylogeny_redo/backtrans_nucleotide_alignments_trimmed/999at7088.faa.reformatted.aln.trimal'  # Replace with your file path
        if remove_requested == False:
            print('[+] Either deletion of codon position which show composition bias despite RY recoding was not enabled or no codon positions show composition bias after recoding. Exiting script.')
            sys.exit()
        elif len(codon_positions_actual)  == 0: # if no codon positions show bias, we don't need to recode
            print('[+] No codon positions show composition bias after RY-recoding. Exiting script.')
            with open(output_file, 'a') as output:
 #           output.write('positions\tchi2_stat\tchr_sq_df\tp_value\n')  # header    
                output.write('Outcome: RY_recoding_successful.\n') 
            sys.exit()
        elif len(codon_positions_actual) == 2: # i.e if two codon positions need to be removed, only one remains - lets get it
                print(f"[+]   Removing codon(s) {codon_positions_actual} from alignment due to base composition bias (p < 0.05)")
                codon_pos_to_keep = set([1,2,3]) - set(codon_positions_actual)
                codon_pos_to_keep = list(codon_pos_to_keep)[0]
                print(f"[+]   Retaining just codon(s) {codon_pos_to_keep}")
                recoded_sequences_with_deletions = keep_every_n_letters(recoded_sequences_dict, codon_pos_to_keep) # keep just the codon position which isnt biased
                with open(output_file, 'a') as output:
                    output.write('Outcome: Removed_codon_positions_which_still_showed_composition_bias_after_RY_recoding.\n')
        elif len(codon_positions_actual) == 3: # i.e. all codon positions show bias, save to warnings file.
                print('[+]   Warning! All three codon positions still show composition bias, even after RY recoding.')
                output_file = output_prefix + '.base_composition_results.tsv'
                with open(output_file, 'a') as output:
                    output.write('Outcome: Base_composition_bias_remains_after_RY_recoding_suggest_discarding_alignment.\n')
                # do something bettre??
        elif set(codon_positions_actual) == set([1]):
                print(f"[+]   Removing codon {codon_positions_actual} from alignment due to base composition bias (p < 0.05)")
                recoded_sequences_with_deletions = interleave_sequences(second_codon_sequences_dict, third_codon_sequences_dict) # order matters, start with first codon positions then second each time
                with open(output_file, 'a') as output:
                    output.write('Outcome: Removed_codon_positions_which_still_showed_composition_bias_after_RY_recoding.\n')
        elif set(codon_positions_actual) == set([2]):
                print(f"[+]   Removing codon {codon_positions_actual} from alignment due to base composition bias (p < 0.05)")
                recoded_sequences_with_deletions = interleave_sequences(first_codon_sequences_dict, third_codon_sequences_dict)
                with open(output_file, 'a') as output:
                    output.write('Outcome: Removed_codon_positions_which_still_showed_composition_bias_after_RY_recoding.\n')
        elif set(codon_positions_actual) == set([3]):
                print(f"[+]   Removing codon {codon_positions_actual} from alignment due to base composition bias (p < 0.05)")
                recoded_sequences_with_deletions = interleave_sequences(second_codon_sequences_dict, second_codon_sequences_dict)
                with open(output_file, 'a') as output:
                    output.write('Outcome: Removed_codon_positions_which_still_showed_composition_bias_after_RY_recoding.\n')
            
        output_alignments_file  = output_prefix + '.recoded.positions.deleted'
        print(f"[+]   Writing RY-recoded and codon position(s) removed sequences to {output_file}")
        if partitions_requested != False:
            parititon_file = output_prefix + '.with_deletions.partitions' 
            print(f"[+]   Writing partition file for the alignment with removed codon position(s) to {parititon_file}")
            number_codons_remaining = 3 - len(codon_positions_actual)
            write_partition_file(recoded_sequences_with_deletions, number_codons_remaining, parititon_file)
        
        with open(output_alignments_file, 'w') as output:
            for seq_id, sequence in recoded_sequences_with_deletions.items():
            # Write each sequence in the desired format
                output.write(f">{seq_id}\n{sequence}\n")
