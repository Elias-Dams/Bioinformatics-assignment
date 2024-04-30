import itertools
import numpy as np
import configparser


def load_config():
    """Load configuration settings from a .ini file."""
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
    """Returns a score based on whether characters 'a' and 'b' match."""
    return match if a == b else mismatch

def initialize_matrix(rows, cols, init_val):
    """Initializes a 2D matrix with specified rows, cols, and initial value."""
    return np.full((rows, cols), init_val)

def pairwise_needleman_wunsch(seq1, seq2, match, mismatch, indel):
    """Performs pairwise alignment using Needleman-Wunsch"""
    rows = len(seq1) + 1
    cols = len(seq2) + 1

    matrix = initialize_matrix(rows, cols, 0)

    # Fill initial row and column for global alignment
    for i in range(1, rows):
        matrix[i, 0] = matrix[i - 1, 0] + indel
    for j in range(1, cols):
        matrix[0, j] = matrix[0, j - 1] + indel

    # Fill remaining matrix
    for i in range(1, rows):
        for j in range(1, cols):
            match_score = matrix[i-1, j-1] + score(seq1[i-1], seq2[j-1], match, mismatch)
            del_score = matrix[i-1, j] + indel
            ins_score = matrix[i, j-1] + indel
            matrix[i, j] = max([match_score, del_score, ins_score])

    rows, cols = matrix.shape
    i, j = rows - 1, cols - 1
    aligned1, aligned2 = "", ""

    while i > 0 or j > 0:
        current = matrix[i, j]
        if i > 0 and j > 0 and current == matrix[i - 1, j - 1] + score(seq1[i - 1], seq2[j - 1], match,
                                                                                 mismatch):
            aligned1 = seq1[i - 1] + aligned1
            aligned2 = seq2[j - 1] + aligned2
            i -= 1
            j -= 1
        elif i > 0 and current == matrix[i - 1, j] + indel:
            aligned1 = seq1[i - 1] + aligned1
            aligned2 = "-" + aligned2
            i -= 1
        elif j > 0 and current == matrix[i, j - 1] + indel:
            aligned1 = "-" + aligned1
            aligned2 = seq2[j - 1] + aligned2
            j -= 1

    return matrix[-1][-1], aligned1, aligned2

def pairwise_smith_waterman(seq1, seq2, match, mismatch, indel):
    """Performs pairwise alignment using Smith Waterman algorithm."""
    rows = len(seq1) + 1
    cols = len(seq2) + 1

    # Initialize the score matrix and the traceback matrix
    score_matrix = initialize_matrix(rows, cols, 0)
    traceback_matrix = initialize_matrix(rows, cols, None)

    # Initialize variables to keep track of the maximum score and its position
    max_score = 0
    max_position = None

    # Fill the score matrix
    for i in range(1, rows):
        for j in range(1, cols):
            match_score = score_matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            del_score = score_matrix[i-1][j] + indel
            ins_score = score_matrix[i][j-1] + indel
            score_matrix[i][j] = max(0, match_score, del_score, ins_score)

            # Update the traceback matrix to keep track of the direction of the maximum score
            if score_matrix[i][j] == 0:
                traceback_matrix[i][j] = None
            elif score_matrix[i][j] == match_score:
                traceback_matrix[i][j] = 'diagonal'
            elif score_matrix[i][j] == del_score:
                traceback_matrix[i][j] = 'up'
            else:
                traceback_matrix[i][j] = 'left'

            # Update the maximum score and its position
            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_position = (i, j)

    # Backtrack to find the alignment
    alignment_seq1 = ''
    alignment_seq2 = ''
    i, j = max_position
    while i > 0 and j > 0 and score_matrix[i][j] > 0:
        if traceback_matrix[i][j] == 'diagonal':
            alignment_seq1 = seq1[i-1] + alignment_seq1
            alignment_seq2 = seq2[j-1] + alignment_seq2
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == 'up':
            alignment_seq1 = seq1[i-1] + alignment_seq1
            alignment_seq2 = '-' + alignment_seq2
            i -= 1
        elif traceback_matrix[i][j] == 'left':
            alignment_seq1 = '-' + alignment_seq1
            alignment_seq2 = seq2[j-1] + alignment_seq2
            j -= 1

    return max_score, alignment_seq1, alignment_seq2

def get_pairwise_alignments(sequences, match, mismatch, indel):
    """Generates pairwise alignments and their scores for all sequence pairs."""
    ids = list(sequences.keys())
    data = list(sequences.values())

    pairwise_combinations1 = itertools.combinations(range(len(data)), 2)
    pairwise_combinations2 = itertools.combinations(range(len(data)), 2)


    print("### Needleman Wunsch:")
    for i, j in pairwise_combinations1:
        seq1, seq2 = data[i], data[j]
        max_score, alignment_seq1, alignment_seq2 = pairwise_needleman_wunsch(seq1, seq2, match, mismatch, indel)
        print(f"Alignment of {ids[i]} vs {ids[j]}:")
        print(f"{ids[i]}: {alignment_seq1}")
        print(f"{ids[j]}: {alignment_seq2}")
        print(f"Score: {max_score}\n")

    print("### Smith Waterman:")
    for i, j in pairwise_combinations2:
        seq1, seq2 = data[i], data[j]
        max_score, alignment_seq1, alignment_seq2 = pairwise_smith_waterman(seq1, seq2, match, mismatch, indel)
        print(f"Alignment of {ids[i]} vs {ids[j]}:")
        print(f"{ids[i]}: {alignment_seq1}")
        print(f"{ids[j]}: {alignment_seq2}")
        print(f"Score: {max_score}\n")


if __name__ == "__main__":
    config = load_config()
    match = int(config['Scores']['match'])
    mismatch = int(config['Scores']['mismatch'])
    indel = int(config['Scores']['indel'])

    sequences = read_fasta('cs_assignment.fasta')
    get_pairwise_alignments(sequences, match=5, mismatch=-2, indel=-4)
