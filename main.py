import configparser
import itertools

def load_config():
    """Load and return the configuration from a .ini file."""
    config = configparser.ConfigParser()
    config.read('config.ini')
    return config

def read_fasta(file_path):
    """Reads sequences from a FASTA file and returns a dictionary mapping sequence IDs to sequences."""
    sequences = {}
    with open(file_path, 'r') as f:
        sequence_id = None
        sequence_data = []

        for line in f:
            line = line.strip()

            if line.startswith('>'):
                if sequence_id and sequence_data:
                    sequences[sequence_id] = ''.join(sequence_data)
                sequence_id = line[1:].strip()
                sequence_data = []
            else:
                sequence_data.append(line)

        if sequence_id and sequence_data:
            sequences[sequence_id] = ''.join(sequence_data)

    return sequences

def score(a, b, match, mismatch):
    """Returns a score based on whether a and b match or not."""
    return match if a == b else mismatch

def initialize_nd_matrix(dimensions, init_val):
    """Recursively initializes an N-dimensional matrix."""
    if len(dimensions) == 1:
        return [init_val for _ in range(dimensions[0])]
    return [initialize_nd_matrix(dimensions[1:], init_val) for _ in range(dimensions[0])]

def needleman_wunsch(sequences, match, mismatch, indel):
    """Performs global alignment using Needleman-Wunsch algorithm on multiple sequences."""
    seq_ids = list(sequences.keys())
    seq_data = [sequences[id] for id in seq_ids]
    num_sequences = len(seq_data)

    lengths = [len(seq) + 1 for seq in seq_data]
    matrix = initialize_nd_matrix(lengths, 0)

    # Fill initial gaps
    for i in range(1, len(seq_data[0]) + 1):
        matrix[i] = matrix[i - 1] + indel

    # Fill the remaining matrix
    for indices in itertools.product(*(range(1, len(seq) + 1) for seq in seq_data)):
        # Determine possible scores for the current indices
        pairwise_combinations = list(itertools.combinations(indices, 2))
        match_score = matrix[tuple(idx - 1 for idx in indices)] + sum(
            score(seq_data[c[0] - 1][indices[c[0]] - 1], seq_data[c[1] - 1][indices[c[1]] - 1], match, mismatch)
            for c in pairwise_combinations
        )
        ins_scores = [matrix[tuple(indices[:i] + (indices[i] - 1,) + indices[i+1:])] + indel for i in range(num_sequences)]
        matrix[indices] = max([match_score] + ins_scores)

    # Backtrack to get alignments
    alignments = [""] * num_sequences
    indices = tuple(len(seq) for seq in seq_data)

    while all(i > 0 for i in indices):
        pairwise_combinations = list(itertools.combinations(indices, 2))
        match_score = matrix[tuple(idx - 1 for idx in indices)] + sum(
            score(seq_data[c[0] - 1][indices[c[0]] - 1], seq_data[c[1] - 1][indices[c[1]] - 1], match, mismatch)
            for c in pairwise_combinations
        )

        ins_scores = [matrix[tuple(indices[:i] + (indices[i] - 1,) + indices[i+1:])] + indel for i in range(num_sequences)]
        max_score = max([match_score] + ins_scores)

        if max_score == match_score:
            for i in range(num_sequences):
                alignments[i] = seq_data[i][indices[i] - 1] + alignments[i]
            indices = tuple(idx - 1 for idx in indices)
        else:
            idx = ins_scores.index(max_score)
            alignments[idx] = "-" + alignments[idx]
            indices = tuple(indices[:idx] + (indices[idx] - 1,) + indices[idx+1:])

    # Output the alignment and its score
    for i in range(num_sequences):
        print(f"{seq_ids[i]}: {alignments[i]}")
    print(f"Alignment Score: {matrix[-1]}")

def smith_waterman(sequences, match, mismatch, indel):
    """Performs local alignment using Smith-Waterman algorithm on multiple sequences."""
    seq_ids = list(sequences.keys())
    seq_data = [sequences[id] for id in seq_ids]
    num_sequences = len(seq_data)

    lengths = [len(seq) + 1 for seq in seq_data]
    matrix = initialize_nd_matrix(lengths, 0)

    # Fill the remaining matrix
    for indices in itertools.product(*(range(1, len(seq) + 1) for seq in seq_data)):
        pairwise_combinations = list(itertools.combinations(indices, 2))
        match_score = matrix[tuple(idx - 1 for idx in indices)] + sum(
            score(seq_data[c[0] - 1][indices[c[0]] - 1], seq_data[c[1] - 1][indices[c[1]] - 1], match, mismatch)
            for c in pairwise_combinations
        )
        ins_scores = [matrix[tuple(indices[:i] + (indices[i] - 1,) + indices[i+1:])] + indel for i in range(num_sequences)]
        matrix[indices] = max([0, match_score] + ins_scores)

    # Find the maximum score and its position
    max_val = 0
    max_pos = (0, 0)

    for indices in itertools.product(*(range(1, len(seq) + 1) for seq in seq_data)):
        if matrix[indices] > max_val:
            max_val = matrix[indices]
            max_pos = indices

    # Backtrack to get the alignment
    indices = max_pos
    alignments = [""] * num_sequences

    while matrix[indices] > 0:
        pairwise_combinations = list(itertools.combinations(indices, 2))
        match_score = matrix[tuple(idx - 1 for idx in indices)] + sum(
            score(seq_data[c[0] - 1][indices[c[0]] - 1], seq_data[c[1] - 1][indices[c[1]] - 1], match, mismatch)
            for c in pairwise_combinations
        )

        ins_scores = [matrix[tuple(indices[:i] + (indices[i] - 1,) + indices[i+1:])] + indel for i in range(num_sequences)]
        max_score = max([match_score] + ins_scores)

        if max_score == match_score:
            for i in range(num_sequences):
                alignments[i] = seq_data[i][indices[i] - 1] + alignments[i]
            indices = tuple(idx - 1 for idx in indices)
        else:
            idx = ins_scores.index(max_score)
            alignments[idx] = "-" + alignments[idx]
            indices = tuple(indices[:idx] + (indices[idx] - 1,) + indices[idx+1:])

    # Output the alignment and its score
    for i in range(num_sequences):
        print(f"{seq_ids[i]}: {alignments[i]}")
    print(f"Alignment Score: {matrix[max_pos]}")

if __name__ == "__main__":
    config = load_config()
    match = int(config['Scores']['match'])
    mismatch = int(config['Scores']['mismatch'])
    indel = int(config['Scores']['indel'])

    sequences = read_fasta('cs_assignment.fasta')
    print("Global Alignment:")
    needleman_wunsch(sequences, match, mismatch, indel)
    print("\nLocal Alignment:")
    smith_waterman(sequences, match, mismatch, indel)
