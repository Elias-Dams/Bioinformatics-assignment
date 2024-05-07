import numpy as np
import itertools


class MultipleSequenceAligner:

    def __init__(self, match: int = 5, mismatch: int = -2, indel: int = -4, two_gaps: int = 0):
        self.match = match
        self.mismatch = mismatch
        self.indel = indel
        self.two_gaps = two_gaps

    def _read_file(self, file_path):
        """Reads sequences from a FASTA file and returns a dictionary mapping sequence IDs to sequences."""
        sequences = {}
        sequence_id = None
        sequence_data = []

        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    if sequence_id is not None:
                        sequences[sequence_id] = ''.join(sequence_data)
                        sequence_data = []
                    sequence_id = line[1:].strip()
                else:
                    sequence_data.append(line.strip())

            if sequence_id is not None:
                sequences[sequence_id] = ''.join(sequence_data)

        return sequences

    def _get_previous_neighbours(self, index: tuple) -> list[tuple[int, ...]]:
        """Generates previous neighbours for a given index, considering all combinations of decrementing each index
        by 0 or 1."""
        num_sequences = len(index)
        offsets = itertools.product([-1, 0], repeat=num_sequences)
        neighbours = [
            neighbour
            for neighbour in (tuple(map(sum, zip(index, offset))) for offset in offsets if any(offset))
            if all(n >= 0 for n in neighbour)
        ]
        return neighbours

    def _get_pairwise_score(self, neighbour: tuple, pair: tuple, sequences: dict, index: tuple[int]) -> int:
        """Calculates the pairwise score between two sequences at given indices based on alignment parameters."""
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

    def _calculate_score_and_direction(self, index: tuple, sequences: dict, score_matrix: np.ndarray, method: str) -> tuple[int, tuple[int] or None]:
        """Calculates the optimal score and the previous cell to move from for a given matrix index."""
        # Initialise the top left as 0
        if all(i == 0 for i in index):
            return 0, None

        # Get all unique pairs of indices for the given dimension
        dimension = len(sequences)
        pairs = [(i, j) for i in range(dimension) for j in range(i + 1, dimension)]

        # Calculate scores for all previous neighbouring cells
        scores = {}
        for neighbour in self._get_previous_neighbours(index):
            total_score = score_matrix[neighbour]
            total_score += sum(self._get_pairwise_score(neighbour, pair, sequences, index) for pair in pairs)
            scores[total_score] = neighbour

        max_score = max(scores)

        if method == "global":
            return max_score, scores[max_score]
        elif method == "local":
            return max(max_score, 0), scores[max_score] if max_score > 0 else None
        else:
            raise ValueError(f"Method must be either 'global' or 'local'. Not '{method}'")

    def _get_current_position(self, sequences: dict, method: str, score_matrix: np.ndarray) -> tuple[int, ...]:
        """Determines the current position in the matrix to start the traceback, based on the alignment method."""
        if method == "local":
            return np.unravel_index(np.argmax(score_matrix), score_matrix.shape)
        return tuple(len(seq) for seq in sequences.values())

    def _update_sequences(self, aligned_sequences: dict, sequences: dict, current_position: tuple, previous_neighbour: tuple) -> None:
        """Updates the aligned sequences based on the current and previous positions in the traceback process."""
        for sequence_index, sequences_id in enumerate(sequences):
            if current_position[sequence_index] == previous_neighbour[sequence_index]:
                aligned_sequences[sequences_id] = "-" + aligned_sequences[sequences_id]
            else:
                char_index = current_position[sequence_index] - 1
                aligned_sequences[sequences_id] = sequences[sequences_id][char_index] + aligned_sequences[sequences_id]

    def _should_break(self, method, previous_neighbour, score_matrix):
        """Determines whether to break out of the alignment loop based on the alignment method and the matrix
        position."""
        if method == "global":
            return all(i == 0 for i in previous_neighbour)
        elif method == "local":
            return np.all(score_matrix[tuple(previous_neighbour)] == 0)
        else:
            raise ValueError(f"Invalid method: {method}")

    def _get_alignment(self, sequences: dict, method: str) -> tuple[dict, float]:
        """Computes the sequence alignment based on the given method, returning the aligned sequences and the total
        score."""
        matrix_dimensions = [len(sequences[id]) + 1 for id in sequences]
        score_matrix = np.zeros(matrix_dimensions)

        prev_cell_dict = {}

        # Populate the matrix with scores.
        # Get all possible indices
        all_indices = list(np.ndindex(score_matrix.shape))
        # Sort indices based on the sum of their elements, so they are sorted diagonally.
        sorted_indices = sorted(all_indices, key=lambda index: sum(index))
        for index in sorted_indices:
            calculated_score, chosen_direction = self._calculate_score_and_direction(index, sequences, score_matrix, method)
            score_matrix[index] = calculated_score
            prev_cell_dict[index] = chosen_direction

        # Backtrack the matrices
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

    def align_sequence(self, input_file: str, method: str = "global", write: bool = False) -> None:
        """Reads sequences from a file and performs alignment using the specified method, then outputs the results."""
        sequences = self._read_file(input_file)

        aligned_sequences, score = self._get_alignment(sequences, method=method)

        # Output the aligned sequences and the score to the console
        output_lines = []
        for key in aligned_sequences:
            line = f"{key}: {aligned_sequences[key]}"
            print(line)
            output_lines.append(line)
        score_line = f"Score: {score}"
        print(score_line)
        output_lines.append(score_line)

        # If write is True, write the output to a text file named based on the input file
        if write:
            output_file = input_file.replace(".fasta","") + "_output.txt"
            with open(output_file, 'w') as f:
                for line in output_lines:
                    f.write(line + "\n")


SequenceAligner = MultipleSequenceAligner(match=5, mismatch=-2, indel=-4, two_gaps=0)
SequenceAligner.align_sequence("cs_assignment.fasta", method="global", write=True)
