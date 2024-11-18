import time

def suffix_array(s):
    """Creates the suffix array for the string."""
    return sorted(range(len(s)), key=lambda i: s[i:])

def bwt_transform(s):
    """Performs the Burrows-Wheeler Transform (BWT) on the string."""
    s += '$'  # Append end-of-string character
    sa = suffix_array(s)
    return ''.join(s[i - 1] for i in sa), sa

def fm_index(bwt):
    """Constructs the FM-index, consisting of the C-array and the occurrence table."""
    c_array = {}
    occ_table = {ch: [0] * (len(bwt) + 1) for ch in "ACGT$"}
    
    # Build the occurrence table
    for i, char in enumerate(bwt):
        for ch in "ACGT$":
            occ_table[ch][i + 1] = occ_table[ch][i] + (1 if char == ch else 0)
    
    # Build the C-array
    total = 0
    for ch in sorted(occ_table.keys()):
        c_array[ch] = total
        total += occ_table[ch][-1]
    
    return c_array, occ_table

def count_occurrences(substring, bwt, c_array, occ_table):
    """Finds all occurrences of the substring using the FM-index and returns the positions."""
    m = len(substring)
    if m == 0:
        return []

    # Start at the last character of the substring
    char = substring[-1]
    if char not in c_array:
        return []

    # Start with the range of the last character
    l = c_array[char]
    r = c_array[char] + occ_table[char][-1]

    # Process the substring backwards
    for i in range(m - 2, -1, -1):
        char = substring[i]
        if char not in c_array:
            return []

        l = c_array[char] + occ_table[char][l]
        r = c_array[char] + occ_table[char][r]

        if l >= r:
            return []

    return list(range(l, r))

def exact_pattern_matching(sequence, pattern):
    """Performs BWT and pattern matching on the given sequence."""
    # Step 1: Perform BWT on the genome sequence
    bwt, sa = bwt_transform(sequence)

    # Step 2: Build the FM-index
    c_array, occ_table = fm_index(bwt)

    # Step 3: Perform exact pattern matching
    positions = count_occurrences(pattern, bwt, c_array, occ_table)

    # Return results
    return len(positions) > 0

def search_in_all_sequences(filename, pattern):
    """Reads the genome sequences from a file, applies BWT, and performs exact pattern matching on each sequence."""
    with open(filename, 'r') as file:
        data = file.read().replace('\n', '').replace(' ', '')  # Clean up the data by removing newlines and spaces

        # Split the genome into sequences (adjust length if needed)
        sequences = [data[i:i + 20] for i in range(0, len(data), 20)]  # Adjust 20 to sequence length

    found_sequences = []
    start_time = time.time()

    # Process each sequence and print C-array, Occ table and perform pattern search
    for seq_num, sequence in enumerate(sequences):
        # Perform BWT on the sequence
        bwt, sa = bwt_transform(sequence)

        # Build FM-index
        c_array, occ_table = fm_index(bwt)

        # Print details for the current sequence
        print(f"\n--- Sequence {seq_num + 1} ---")
        print(f"Original Sequence: {sequence}")
        print(f"BWT: {bwt}")
        print(f"C-array: {c_array}")
        print(f"Occurrence Table: {occ_table}")

        # Search for the pattern in the current sequence
        if exact_pattern_matching(sequence, pattern):
            found_sequences.append(seq_num + 1)  # Sequence number is 1-based

    end_time = time.time()

    # Calculate efficiency
    search_time = end_time - start_time

    # Output results
    return found_sequences, len(found_sequences), search_time

if __name__ == "_main_":
    filename = 'corona.txt'  # Specify your input file here
    pattern = "AGGCTA"  # This is the pattern we want to find in each sequence

    # Perform pattern matching on all sequences
    found_sequences, total_found, search_time = search_in_all_sequences(filename, pattern)

    # Output the results
    print(f"\n--- Summary ---")
    print(f"Total number of sequences where pattern '{pattern}' was found: {total_found}")
    if found_sequences:
        print(f"Pattern '{pattern}' found in sequence numbers: {found_sequences}")
    else:
        print(f"Pattern '{pattern}' not found in any sequence.")
    print(f"Time taken for pattern search: {search_time:.4f} seconds")