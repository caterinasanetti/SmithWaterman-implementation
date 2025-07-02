# SmithWaterman-implementation
This code implements the Smith Waterman algorithm , with the possibility of using also affine gap penalities 
# Python Local Sequence Aligner

A command-line tool for performing local sequence alignment using dynamic programming (Smith-Waterman algorithm) with support for both linear and affine gap penalties.

## Features

*   **Local Alignment:** Implements the Smith-Waterman algorithm to find the most similar regions between two sequences.
*   **Gap Penalty Modes:** Supports both linear and affine gap penalties.
*   **Custom Scoring:** Allows specifying scores for matches, mismatches, and gap penalties.
*   **Multiple Alignments:** Identifies and can display multiple high-scoring local alignments, starting from the highest scoring ones.
*   **Detailed Output:** Provides aligned sequences, score, start positions, length, and alignment statistics (matches, mismatches, gaps).
*   **Visualization:** Includes a simple visualization line indicating matches, mismatches, and gaps.
*   **Command-Line Interface:** Easy to use via the terminal.

## How it Works

The tool uses dynamic programming to fill a scoring matrix based on the provided sequences and scoring parameters. The core algorithm is a variation of the Smith-Waterman algorithm for local alignment, which involves:

1.  Initializing a score matrix where the first row and column are filled with zeros.
2.  Iterating through the matrix, calculating the score for each cell based on the scores of adjacent cells and the match/mismatch score for the corresponding characters. The score is the maximum of:
    *   Zero (enabling local alignment - any negative score becomes zero, allowing new alignments to start).
    *   The score from the diagonal cell plus the match/mismatch score.
    *   The score from the cell above plus the gap penalty (or gap extension penalty).
    *   The score from the cell to the left plus the gap penalty (or gap extension penalty).
3.  Storing the decision made at each cell (which neighbor resulted in the maximum score) in a traceback matrix.
4.  Once the matrix is filled, finding the cell(s) with the highest score.
5.  Performing traceback from the highest scoring cell(s) (or all positive-scoring cells, sorted by score) back through the matrix following the recorded decisions until a cell with a score of zero is reached. This process reconstructs the alignment.
6.  For **affine gap penalties**, the implementation uses additional matrices (similar to Gotoh's algorithm) to track the cost of opening and extending gaps separately.

## Installation

This script requires Python 3. It uses only standard Python libraries (`argparse`, `collections`).

1.  Save the provided code as a Python file, for example, `local_aligner.py`.
2.  Make sure you have Python 3 installed.

## Usage
After cloning the reopsitory to your local machine , run the script from your terminal, providing the two sequences as command-line arguments.


`python3 local_aligner.py <sequence1> <sequence2> [options]`


