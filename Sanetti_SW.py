import argparse #Used for parsing command-line arguments.
from collections import deque # Used for efficient traceback path exploration.


class LocalAligner:
    """Performs local sequence alignment using either linear or affine gap penalties."""

    def __init__(self, seq1, seq2, match=1, mismatch=-1, gap=-2, mode='linear', gap_open=None, gap_extend=-1):
        """
        Initialize aligner with sequences and scoring parameters.

        Args:
            seq1, seq2: Input sequences to align
            match: Score for matching characters
            mismatch: Penalty for mismatches
            gap: Linear gap penalty or affine gap open penalty
            mode: 'linear' or 'affine' gap penalty mode
            gap_open: Affine gap open penalty (overrides gap if provided)
            gap_extend: Affine gap extension penalty
        """

        self.seq1 = seq1
        self.seq2 = seq2
        self.match = match
        self.mismatch = mismatch
        self.mode = mode

        # Set gap penalties based on mode
        if mode == 'affine':
            self.gap_open = gap_open if gap_open is not None else gap
            self.gap_extend = gap_extend
        else:
            self.gap_penalty = gap  # Simple linear gap penalty

        # Initialize score matrices
        self._init_matrices()

    def _init_matrices(self):
        """Initialize the scoring and traceback matrices."""
        rows, cols = len(self.seq1) + 1, len(self.seq2) + 1

        # Main score matrix stores alignment scores
        self.score = [[0] * cols for _ in range(rows)]

        # For affine gaps: matrices for gaps in seq1 (P) and seq2 (Q)
        if self.mode == 'affine':
            self.P = [[0] * cols for _ in range(rows)]  # Vertical gaps
            self.Q = [[0] * cols for _ in range(rows)]  # Horizontal gaps

        # Traceback matrix stores directions (diagonal, up, left)
        self.traceback = [[[] for _ in range(cols)] for _ in range(rows)]

        self.max_score = 0  #Tracks the highest alignment score
        self.max_positions = [] # Positions of maximum scores
        self.all_positions = []  # All positions with score > 0

    def _match_score(self, a, b):
        """Return score for aligning two characters."""
        return self.match if a == b else self.mismatch

    def align(self):
        """Fill the scoring matrices using dynamic programming.
        Steps :
        1) Iterates over each cell (i, j) in the Main Score matrix.
        2)Calls _fill_cell(i, j) to compute the score
        3)Tracks the maximum score and positions where it occurs.
        4)Collects all positions with a score > 0.
        5) Sorts all positions by score in descending order."""


        for i in range(1, len(self.seq1) + 1):
            for j in range(1, len(self.seq2) + 1):
                # Calculate score for current cell
                self._fill_cell(i, j)

                # Track maximum scores
                if self.score[i][j] > self.max_score:
                    self.max_score = self.score[i][j]
                    self.max_positions = [(i, j)]
                elif self.score[i][j] == self.max_score:
                    self.max_positions.append((i, j))

                # Track all positions with positive score
                if self.score[i][j] > 0:
                    self.all_positions.append((i, j))

        # Sort all positions by score (highest first)
        self.all_positions.sort(key=lambda p: -self.score[p[0]][p[1]])

    def _fill_cell(self, i, j):
        """Calculate score for cell (i,j).
        For affine gaps :
        - Uses separate matrices P and Q for gap penalties
        - Computes P[i][j] (gap in seq1) , Q[i][j] (gap in seq2) and score[i][j] = max(0, diag, P[i][j], Q[i][j])

        For linear gaps :
        - Computes score[i][j] = max(0, diag, up + gap , left + gap)
        """
        # Score for matching current characters
        match = self._match_score(self.seq1[i - 1], self.seq2[j - 1])
        diag = self.score[i - 1][j - 1] + match

        if self.mode == 'affine':
            # Affine gap penalties
            self.P[i][j] = max(
                self.score[i - 1][j] + self.gap_open + self.gap_extend,
                self.P[i - 1][j] + self.gap_extend
            )
            self.Q[i][j] = max(
                self.score[i][j - 1] + self.gap_open + self.gap_extend,
                self.Q[i][j - 1] + self.gap_extend
            )
            self.score[i][j] = max(0, diag, self.P[i][j], self.Q[i][j])
        else:
            # linear gap penalty
            up = self.score[i - 1][j] + self.gap_penalty
            left = self.score[i][j - 1] + self.gap_penalty
            self.score[i][j] = max(0, diag, up, left)

        # Record traceback directions
        self._set_traceback(i, j, diag)

    def _set_traceback(self, i, j, diag):
        """Record traceback directions."""
        directions = []
        if self.score[i][j] == 0:
            directions.append(0)  # Stop traceback
        if self.score[i][j] == diag:
            directions.append(1)  # Diagonal (match/mismatch)

        if self.mode == 'affine':
            if self.score[i][j] == self.P[i][j]:
                directions.append(2)  # Up (gap in seq1)
            if self.score[i][j] == self.Q[i][j]:
                directions.append(3)  # Left (gap in seq2)
        else:
            if self.score[i][j] == self.score[i - 1][j] + self.gap_penalty:
                directions.append(2)  # Up
            if self.score[i][j] == self.score[i][j - 1] + self.gap_penalty:
                directions.append(3)  # Left

        self.traceback[i][j] = directions


class AlignmentTraceback:
    """Handles traceback operations to find optimal alignments."""

    def __init__(self, score_matrix, traceback_matrix, seq1, seq2):
        """Initializes with score matrix, traceback matrix, and sequences."""
        self.score = score_matrix
        self.traceback = traceback_matrix
        self.seq1 = seq1
        self.seq2 = seq2

    def find_alignments(self, start_pos):
        """Find all  optimal alignments .
        Steps :
        1)Starts from (i, j) (given position).
        2)Uses a queue to explore all traceback paths.
        3)Reconstructs alignments in reverse order."""
        queue = deque([(start_pos[0], start_pos[1], [], [])])
        alignments = []

        while queue:
            i, j, aln1, aln2 = queue.popleft()

            # Reached beginning of alignment
            if self.score[i][j] == 0:
                alignment = self._create_alignment(aln1, aln2, start_pos)
                alignments.append(alignment)
                continue

            # Explore all possible paths
            for direction in self.traceback[i][j]:
                new_aln1, new_aln2 = aln1.copy(), aln2.copy()

                if direction == 1:  # Diagonal
                    new_aln1.append(self.seq1[i - 1])
                    new_aln2.append(self.seq2[j - 1])
                    queue.append((i - 1, j - 1, new_aln1, new_aln2))
                elif direction == 2:  # Up (gap in seq2)
                    new_aln1.append(self.seq1[i - 1])
                    new_aln2.append('-')
                    queue.append((i - 1, j, new_aln1, new_aln2))
                elif direction == 3:  # Left (gap in seq1)
                    new_aln1.append('-')
                    new_aln2.append(self.seq2[j - 1])
                    queue.append((i, j - 1, new_aln1, new_aln2))

        return alignments

    def _create_alignment(self, aln1, aln2, start_pos):
        """Formats the alignment result and calculate alignment statistics."""
        # Reverse the aligned sequences (traceback builds them backwards)
        aln1 = ''.join(reversed(aln1))
        aln2 = ''.join(reversed(aln2))

        # Calculate alignment statistics
        matches = sum(1 for a, b in zip(aln1, aln2) if a == b and a != '-')
        mismatches = sum(1 for a, b in zip(aln1, aln2) if a != b and a != '-' and b != '-')
        gaps = sum(1 for a, b in zip(aln1, aln2) if a == '-' or b == '-')

        return {
            'seq1': aln1,
            'seq2': aln2,
            'score': self.score[start_pos[0]][start_pos[1]],
            'start_i': start_pos[0],
            'start_j': start_pos[1],
            'matches': matches,
            'mismatches': mismatches,
            'gaps': gaps
        }


def print_alignment(alignment, num):
    """Print formatted alignment output."""
    # Generate visualization line showing matches/mismatches
    vis = []
    for a, b in zip(alignment['seq1'], alignment['seq2']):
        if a == b and a != '-':
            vis.append('|')  # Match
        elif a == '-' or b == '-':
            vis.append(' ')  # Gap
        else:
            vis.append('*')  # Mismatch

    alignment_length = len(alignment['seq1'])
    print(f"\nAlignment {num}:")
    print(f"Start position: seq1={alignment['start_i']}, seq2={alignment['start_j']}")
    print(f"Score: {alignment['score']}")
    print(f"Length: {alignment_length}")
    print(alignment['seq1'])
    print(''.join(vis))
    print(alignment['seq2'])
    print(f" {alignment['matches']} matches, {alignment['mismatches']} mismatches, {alignment['gaps']} gaps")


def main():
    """Command-line interface for the sequence aligner."""
    parser = argparse.ArgumentParser(description="Local sequence alignment with linear/affine gap penalties")
    parser.add_argument("seq1", help="First sequence")
    parser.add_argument("seq2", help="Second sequence")
    parser.add_argument("--match", type=int, default=1, help="Match score (default: 1)")
    parser.add_argument("--mismatch", type=int, default=-1, help="Mismatch penalty (default: -1)")
    parser.add_argument("--gap", type=int, default=-2, help="Gap penalty (linear) or gap open (affine) (default: -2)")
    parser.add_argument("--gap_open", type=int, help="Affine gap open penalty (overrides --gap)")
    parser.add_argument("--gap_extend", type=int, default=-1, help="Affine gap extension penalty (default: -1)")
    parser.add_argument("--mode", choices=['linear', 'affine'], default='linear',
                        help="Alignment mode (default: linear)")
    parser.add_argument("--limit", type=int, help="Maximum number of alignments to show")
    args = parser.parse_args()

    # Initialize aligner with command-line arguments
    aligner = LocalAligner(
        args.seq1, args.seq2,
        match=args.match,
        mismatch=args.mismatch,
        gap=args.gap,
        mode=args.mode,
        gap_open=args.gap_open,
        gap_extend=args.gap_extend
    )

    # Perform alignment
    aligner.align()

    # Initialize traceback handler
    traceback = AlignmentTraceback(
        aligner.score,
        aligner.traceback,
        args.seq1,
        args.seq2
    )

    # Print summary information
    print(f"\nAlignment mode: {args.mode}")
    print(f"Optimal score: {aligner.max_score}")
    print(f"Found {len(aligner.all_positions)} alignments with score > 0")

    # Process and display alignments
    seen_alignments = set()
    alignment_count = 0

    for pos in aligner.all_positions:
        for alignment in traceback.find_alignments(pos):
            # Ensure each alignment is only shown once
            alignment_key = (alignment['start_i'], alignment['start_j'],
                             alignment['seq1'], alignment['seq2'])

            if alignment_key not in seen_alignments:
                seen_alignments.add(alignment_key)
                alignment_count += 1
                print_alignment(alignment, alignment_count)

                # Stop if we've reached the limit
                if args.limit and alignment_count >= args.limit:
                    return


if __name__ == "__main__":
    main()