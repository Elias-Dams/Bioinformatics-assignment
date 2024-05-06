from typing import List, Tuple

import numpy as np
import itertools


class MultiSequenceAligner:

    def __init__(self, match: int = 5, mismatch: int = -2, indel: int = -4, two_gaps: int = 0):
        self.match = match
        self.mismatch = mismatch
        self.indel = indel
        self.two_gaps = two_gaps

    def _read_file(self, file_path):
        """Reads sequences from a FASTA file and returns a dictionary mapping sequence IDs to sequences."""
        sequences = {}
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    if 'sequence_data' in locals():
                        sequences[sequence_id] = ''.join(sequence_data)
                    sequence_id = line[1:].strip()
                    sequence_data = []
                else:
                    sequence_data.append(line.strip())
            if 'sequence_data' in locals():
                sequences[sequence_id] = ''.join(sequence_data)
        return sequences

    def _get_previous_neighbours(self, index: tuple[int]) -> list[tuple[int, ...]]:
        num_sequences = len(index)
        # Generate potential offset combinations without the all-zero case
        offsets = itertools.product([-1, 0], repeat=num_sequences)
        neighbours = [tuple(map(sum, zip(index, offset))) for offset in offsets if any(offset)]
        return neighbours

    def _get_pairwise_score(self, neighbour: tuple[int, ...], pair: tuple[int, int], sequences: dict, index: tuple[int]) -> int:
        i, j = pair
        neighbour_i, neighbour_j = neighbour[i], neighbour[j]
        idx_i, idx_j = index[i], index[j]
        diff_i, diff_j = idx_i - neighbour_i, idx_j - neighbour_j
        id_i, id_j = list(sequences.keys())[i], list(sequences.keys())[j]

        # handle the 4 cases of pairwise alignment
        if diff_i == 1 and diff_j == 1:
            aa_i = sequences[id_i][idx_i - 1]
            aa_j = sequences[id_j][idx_j - 1]
            score = self.match if aa_i == aa_j else self.mismatch
        elif diff_i == 1 and diff_j == 0:
            score = self.indel
        elif diff_i == 0 and diff_j == 1:
            score = self.indel
        elif diff_i == 0 and diff_j == 0:
            score = self.two_gaps
        else:
            raise ValueError(f"Invalid offset: {diff_i}, {diff_j}")

        return score

    def _calculate_score_and_direction(self, index: tuple[int], sequences: dict, matrix: np.ndarray, method: str) -> tuple[int, tuple[int] or None]:
        """Calculates the score and direction for an index in de matrix given the sequences."""
        # set the initial values of the matrix
        if any(i == 0 for i in index):
            if method == "local":
                return 0, None
            elif method == "global":
                return self.indel * max(index), None
            else:
                raise ValueError(f"Method must be either 'global' or 'local'. Not '{method}'")

        # Get all unique pairs of indices for the given dimension
        dimension = len(sequences)
        pairs = [(i, j) for i in range(dimension) for j in range(i + 1, dimension)]

        # Calculate scores for all previous neighbouring cells
        scores = {}
        for neighbour in self._get_previous_neighbours(index):
            total_score = matrix[neighbour]
            total_score += sum(self._get_pairwise_score(neighbour, pair, sequences, index) for pair in pairs)
            scores[total_score] = neighbour

        max_score = max(scores)

        if method == "global":
            return max_score, scores[max_score]
        elif method == "local":
            return max(max_score, 0), scores[max_score] if max_score > 0 else None
        else:
            raise ValueError(f"Method must be either 'global' or 'local'. Not '{method}'")

    def _get_alignment(self, sequences: dict, method: str) -> list[str]:
        matrix_dimensions = [len(sequences[id]) + 1 for id in sequences]
        score_matrix = np.zeros(matrix_dimensions)

        prev_cell_dict = {}

        # Populate the matrices with scores.
        # Get all possible indices
        all_indices = list(np.ndindex(score_matrix.shape))
        # Sort indices based on the sum of their elements, so they are sorted diagonally.
        sorted_indices = sorted(all_indices, key=lambda index: sum(index))
        for index in sorted_indices:
            calculated_score, chosen_direction = self._calculate_score_and_direction(index, sequences, score_matrix,
                                                                                     method)
            score_matrix[index] = calculated_score
            prev_cell_dict[index] = chosen_direction

        # score_matrix, prev_cell_dict
        # ....
        return ["", ""]

    def align_sequence(self, input_file, method="global"):
        sequences = self._read_file(input_file)

        aligned_sequences = self._get_alignment(sequences, method=method)


SequenceAligner = MultiSequenceAligner(match=5, mismatch=-2, indel=-4, two_gaps=0)
SequenceAligner.align_sequence("cs_assignment.fasta", method="global")
