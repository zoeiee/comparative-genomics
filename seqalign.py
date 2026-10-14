#!/usr/bin/env python3
"""Needleman-Wunsch global alignment and edit-distance tools for PS1 Q1."""
import argparse, csv
from pathlib import Path


def read_fasta(path):
    """Read one or more FASTA records and concatenate their sequence lines."""
    seq = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(('>', ';', '#')):
                continue
            seq.append(''.join(line.split()).upper())
    result = ''.join(seq)
    if not result:
        raise ValueError(f"No sequence found in FASTA/text file: {path}")
    return result


def needleman_wunsch(a, b, match=1, mismatch=-1, gap=-1):
    """Return (optimal_score, aligned_a, aligned_b, score_matrix, traceback_matrix).

    Ties are resolved diagonal, then up, then left. Any tied path is optimal.
    """
    n, m = len(a), len(b)
    F = [[0] * (m + 1) for _ in range(n + 1)]
    T = [[None] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        F[i][0], T[i][0] = i * gap, 'U'
    for j in range(1, m + 1):
        F[0][j], T[0][j] = j * gap, 'L'
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            sub = match if a[i - 1] == b[j - 1] else mismatch
            choices = (F[i-1][j-1] + sub, F[i-1][j] + gap, F[i][j-1] + gap)
            best = max(choices)
            F[i][j] = best
            T[i][j] = 'D' if choices[0] == best else ('U' if choices[1] == best else 'L')
    i, j = n, m
    aa, bb = [], []
    while i or j:
        move = T[i][j]
        if move == 'D':
            aa.append(a[i-1]); bb.append(b[j-1]); i -= 1; j -= 1
        elif move == 'U':
            aa.append(a[i-1]); bb.append('-'); i -= 1
        elif move == 'L':
            aa.append('-'); bb.append(b[j-1]); j -= 1
        else:
            raise RuntimeError(f"Invalid traceback at ({i},{j})")
    return F[n][m], ''.join(reversed(aa)), ''.join(reversed(bb)), F, T


def edit_distance(a, b):
    """Levenshtein distance: min substitutions + insertions + deletions."""
    # Keep only two rows to reduce memory to O(min(len(a), len(b))).
    if len(b) > len(a):
        a, b = b, a
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i] + [0] * len(b)
        for j, cb in enumerate(b, 1):
            cur[j] = min(prev[j] + 1, cur[j-1] + 1, prev[j-1] + (ca != cb))
        prev = cur
    return prev[-1]


def write_matrices(path, a, b, F, T):
    """Write score and traceback matrices into one readable CSV file."""
    path = Path(path)
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['SCORE MATRIX'] + ['-'] + list(b))
        for i, row in enumerate(F):
            label = '-' if i == 0 else a[i-1]
            w.writerow([label] + row)
        w.writerow([])
        w.writerow(['TRACEBACK MATRIX'] + ['-'] + list(b))
        for i, row in enumerate(T):
            label = '-' if i == 0 else a[i-1]
            w.writerow([label] + [x or '' for x in row])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('seq1', help='FASTA file for sequence 1')
    p.add_argument('seq2', help='FASTA file for sequence 2')
    p.add_argument('--match', type=int, default=1, help='match score (default: 1; set to course skeleton value)')
    p.add_argument('--mismatch', type=int, default=-1, help='mismatch score (default: -1; set to course skeleton value)')
    p.add_argument('--gap', type=int, default=-1, help='linear gap score (default: -1; set to course skeleton value)')
    p.add_argument('--matrix-out', help='write score/traceback matrices to CSV')
    p.add_argument('--distance', action='store_true', help='also compute nonnegative Levenshtein distance')
    args = p.parse_args()
    a, b = read_fasta(args.seq1), read_fasta(args.seq2)
    score, aa, bb, F, T = needleman_wunsch(a, b, args.match, args.mismatch, args.gap)
    print('Alignment 1:', aa)
    print('Alignment 2:', bb)
    print('Optimal global alignment score:', score)
    if args.distance:
        print('Levenshtein edit distance:', edit_distance(a, b))
    if args.matrix_out:
        write_matrices(args.matrix_out, a, b, F, T)
        print('Matrices written to:', args.matrix_out)

if __name__ == '__main__':
    main()
