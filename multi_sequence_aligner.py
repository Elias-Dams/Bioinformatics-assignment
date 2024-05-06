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
        offsets = itertools.product([-1, 0], repeat=num_sequences)
        neighbours = [
            neighbour
            for neighbour in (tuple(map(sum, zip(index, offset))) for offset in offsets if any(offset))
            if all(n >= 0 for n in neighbour)
        ]
        return neighbours

    def _get_pairwise_score(self, neighbour: tuple[int, ...], pair: tuple[int, int], sequences: dict, index: tuple[int]) -> int:
        i, j = pair
        neighbour_i, neighbour_j = neighbour[i], neighbour[j]
        index_i, index_j = index[i], index[j]
        diff_i, diff_j = index_i - neighbour_i, index_j - neighbour_j
        id_i, id_j = list(sequences.keys())[i], list(sequences.keys())[j]

        # handle the 4 cases of pairwise alignment
        if diff_i == 1 and diff_j == 1:
            aa_i = sequences[id_i][index_i - 1]
            aa_j = sequences[id_j][index_j - 1]
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
        if all(i == 0 for i in index):
            return 0, None

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

    def _get_current_position(self, sequences, method, score_matrix):
        if method == "local":
            return np.unravel_index(np.argmax(score_matrix), score_matrix.shape)
        return tuple(len(seq) for seq in sequences.values())

    def _update_sequences(self, aligned_sequences, sequences, current_position, previous_neighbour):
        for sequence_index, sequences_id in enumerate(sequences):
            if current_position[sequence_index] == previous_neighbour[sequence_index]:
                aligned_sequences[sequences_id] = "-" + aligned_sequences[sequences_id]
            else:
                char_index = current_position[sequence_index] - 1
                aligned_sequences[sequences_id] = sequences[sequences_id][char_index] + aligned_sequences[sequences_id]

    def _should_break(self, method, previous_neighbour, score_matrix):
        if method == "global":
            return all(i == 0 for i in previous_neighbour)
        elif method == "local":
            return np.all(score_matrix[tuple(previous_neighbour)] == 0)
        else:
            raise ValueError(f"Invalid method: {method}")

    def _get_alignment(self, sequences: dict, method: str) -> tuple[dict, float]:
        matrix_dimensions = [len(sequences[id]) + 1 for id in sequences]
        score_matrix = np.zeros(matrix_dimensions)

        prev_cell_dict = {}

        # Populate the matrices with scores.
        # Get all possible indices
        all_indices = list(np.ndindex(score_matrix.shape))
        # Sort indices based on the sum of their elements, so they are sorted diagonally.
        sorted_indices = sorted(all_indices, key=lambda index: sum(index))
        for index in sorted_indices:
            calculated_score, chosen_direction = self._calculate_score_and_direction(index, sequences, score_matrix, method)
            score_matrix[index] = calculated_score
            prev_cell_dict[index] = chosen_direction

        aligned_sequences = {key: "" for key in sequences}
        current_position = self._get_current_position(sequences, method, score_matrix)
        total_score = float(score_matrix[current_position])

        while True:
            previous_neighbour = prev_cell_dict[current_position]
            self._update_sequences(aligned_sequences, sequences, current_position, previous_neighbour)

            if self._should_break(method, previous_neighbour, score_matrix):
                break

            current_position = previous_neighbour

        return aligned_sequences, total_score

    def align_sequence(self, input_file, method="global"):
        sequences = self._read_file(input_file)

        aligned_sequences, score = self._get_alignment(sequences, method=method)

        for key in aligned_sequences:
            print(f"{key}: {aligned_sequences[key]}")
        print(f"Score: {score}")


SequenceAligner = MultiSequenceAligner(match=5, mismatch=-2, indel=-4, two_gaps=0)
SequenceAligner.align_sequence("cs_assignment.fasta", method="global")
