import re


def find_target_sequences(rna_sequence, pattern):
    target_sequences = []
    #define pattern
    converted_pattern = pattern.replace('N', '[ACGU]').replace('V', '[ACG]')
    pattern = re.compile(converted_pattern)
    
    for i in range(len(rna_sequence) - 28):
        candidate = rna_sequence[i:i+29]
        if pattern.fullmatch(candidate[-5:]):
            target_sequences.append(candidate)
    
    return target_sequences

def reverse_complement(dna_sequence):
    # Define the complement rules for DNA nucleotides
    complement = {
        'A': 'U',
        'U': 'A',
        'C': 'G',
        'G': 'C'
    }

    # Reverse the DNA sequence
    reversed_sequence = dna_sequence[::-1]

    # Get the complement of the reversed sequence
    reverse_complement_sequence = ''.join(complement[base] for base in reversed_sequence)

    return reverse_complement_sequence

def main():
    # Get the RNA sequence
    rna_sequence = input("Enter the RNA sequence (5' => 3'): ").strip().upper()
    pattern = input("Enter the pattern: ").strip().upper()
    # Validate RNA sequence
    if not re.fullmatch(r'[ACGU]+', rna_sequence):
        print("Invalid RNA sequence. Please enter a sequence containing only A, C, G, U.")
        return
    
    target_sequences = find_target_sequences(rna_sequence, pattern)
    #get the final gRNA
    gRNA_seq = [reverse_complement(seq[:-5]) for seq in target_sequences]



    if target_sequences:
        print("Found target sequences:")
        for seq in gRNA_seq:
            print(seq)
    else:
        print("No target sequences found.")

if __name__ == "__main__":
    main()
