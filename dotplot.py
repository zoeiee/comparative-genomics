#!/usr/bin/env python3
"""K-mer/spaced-seed dotplots and inversion-aware matching for PS1 Q2."""
import argparse, itertools, json
from pathlib import Path
from seqalign import read_fasta


def index_seeds(seq, span, step=1, start_limit=None):
    """Map each sampled seed key to all start positions in seq."""
    last = len(seq) - span + 1
    if start_limit is not None:
        last = min(last, start_limit)
    idx = {}
    offsets = tuple(range(0, span, step))
    for i in range(max(0, last)):
        key = ''.join(seq[i + off] for off in offsets)
        idx.setdefault(key, []).append(i)
    return idx, offsets


def find_seed_hits(seq1, seq2, span=30, step=1):
    """Return (x,y) matches where sampled bases at offsets 0,step,2*step,... agree."""
    idx, offsets = index_seeds(seq1, span, step)
    hits = []
    for y in range(len(seq2) - span + 1):
        key = ''.join(seq2[y + off] for off in offsets)
        for x in idx.get(key, ()):
            hits.append((x, y))
    return hits


def reverse_complement(seq):
    return seq.translate(str.maketrans('ACGTNacgtn', 'TGCANtgcan'))[::-1]


def find_inversion_hits(seq1, seq2, span=30, step=1):
    """Find matches to reverse-complement(seq1), mapping x back to original seq1.

    Returned x coordinates are in the original sequence. Such hits trace negative-slope
    diagonals when the region is inverted relative to seq2.
    """
    rc = reverse_complement(seq1)
    rev_hits = find_seed_hits(rc, seq2, span, step)
    return [(len(seq1) - (x_rc + span), y) for x_rc, y in rev_hits]


def find_100mers_local_mismatch_rule(seq1, seq2):
    """Find 100-mer pairs with <=2 mismatches in EVERY contiguous 6-base block.

    Candidate generation: the first six bases of every valid pair must contain at least
    four matching positions. Index all 15 choices of four positions from a 6-base window,
    then verify candidates using the full 100-base rule. This is a seed-and-verify method.
    It can still produce many candidates in repetitive genomes, so use only when needed.
    """
    span = 100
    if len(seq1) < span or len(seq2) < span:
        return []
    masks = list(itertools.combinations(range(6), 4))
    candidates = set()
    # Process one mask at a time to avoid retaining 15 genome-sized indexes.
    for mask in masks:
        idx = {}
        for i in range(len(seq1) - span + 1):
            key = ''.join(seq1[i+p] for p in mask)
            idx.setdefault(key, []).append(i)
        for j in range(len(seq2) - span + 1):
            key = ''.join(seq2[j+p] for p in mask)
            for i in idx.get(key, ()):
                candidates.add((i, j))
    hits = []
    for i, j in candidates:
        good = True
        for off in range(span - 5):
            mismatches = sum(seq1[i+off+k] != seq2[j+off+k] for k in range(6))
            if mismatches > 2:
                good = False
                break
        if good:
            hits.append((i, j))
    return sorted(hits)


def plot_hits(hits, n1, n2, output, title, reverse_hits=None):
    try:
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise SystemExit('Plotting needs matplotlib. Install with: python -m pip install matplotlib') from e
    fig, ax = plt.subplots(figsize=(8, 7))
    if hits:
        ax.scatter([h[0] for h in hits], [h[1] for h in hits], s=1, alpha=.45, label='Forward matches')
    if reverse_hits:
        ax.scatter([h[0] for h in reverse_hits], [h[1] for h in reverse_hits], s=2, alpha=.65, label='Reverse-complement matches')
        ax.legend(loc='best')
    ax.set(xlim=(0, n1), ylim=(0, n2), xlabel='Sequence 1 coordinate', ylabel='Sequence 2 coordinate', title=title)
    ax.set_aspect('auto')
    fig.tight_layout()
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('seq1', help='first FASTA file (human in the assignment)')
    p.add_argument('seq2', help='second FASTA file (mouse in the assignment)')
    p.add_argument('plot', help='output plot filename, e.g. results/dotplot.png')
    p.add_argument('--span', type=int, default=30, help='window span: 30, 100, 60, 90, or 120')
    p.add_argument('--step', type=int, default=1, help='sample every step-th base; use 2/3/4 for spaced seeds')
    p.add_argument('--diag-tolerance', type=int, default=None, help='absolute |x-y| tolerance for near-diagonal hits; default 1%% of longer sequence')
    p.add_argument('--inversions', action='store_true', help='also search reverse-complement orientation')
    p.add_argument('--summary-out', help='write hit counts and specificity summary as JSON')
    p.add_argument('--local-mismatch-rule', action='store_true', help='run the Q2(b)(v) 100-mer <=2 mismatches per every 6-base block test')
    args = p.parse_args()
    a, b = read_fasta(args.seq1), read_fasta(args.seq2)
    if args.local_mismatch_rule:
        hits = find_100mers_local_mismatch_rule(a, b)
        span, step = 100, 1
    else:
        if args.span < 1 or args.step < 1 or args.step > args.span:
            p.error('--span and --step must be positive, with step <= span')
        hits = find_seed_hits(a, b, args.span, args.step)
        span, step = args.span, args.step
    reverse_hits = find_inversion_hits(a, b, span, step) if args.inversions else []
    tol = args.diag_tolerance if args.diag_tolerance is not None else max(1, int(.01 * max(len(a), len(b))))
    near = sum(abs(x-y) <= tol for x, y in hits)
    summary = {
        'sequence1_length': len(a), 'sequence2_length': len(b), 'span': span,
        'sample_step': step, 'sampled_bases': len(range(0, span, step)),
        'forward_hits': len(hits), 'near_diagonal_tolerance': tol,
        'near_diagonal_hits': near, 'near_diagonal_percent': (100*near/len(hits) if hits else 0.0),
        'reverse_complement_hits': len(reverse_hits),
        'note': 'Near-diagonal threshold is configurable; compare runs using the same threshold.'
    }
    Path(args.plot).parent.mkdir(parents=True, exist_ok=True)
    plot_hits(hits, len(a), len(b), args.plot, f'span={span}, step={step}', reverse_hits)
    if args.summary_out:
        Path(args.summary_out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.summary_out).write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps(summary, indent=2))
    print('Plot written to:', args.plot)

if __name__ == '__main__':
    main()
