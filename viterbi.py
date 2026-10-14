#!/usr/bin/env python3
"""Two-state GC-content HMM and Viterbi decoder for PS1 Q3."""
import argparse, json, math
from pathlib import Path
from collections import Counter

STATES = ('H', 'L')  # H = high-GC; L = low-GC
BASES = 'ACGT'


def read_hmm_data(path):
    """Read a DNA sequence and +/− annotation from a plain text data file.

    Accepts separate sequence/annotation lines, FASTA-like sequence lines, or whitespace-
    separated tokens. '+' denotes High-GC and '-' denotes Low-GC, per the assignment.
    Header/comment lines are ignored. Returns (sequence, annotation_or_None).
    """
    seq_parts, ann_parts = [], []
    for raw in Path(path).read_text(encoding='utf-8').splitlines():
        line = ''.join(raw.split())
        if not line or line.startswith(('>', '#', ';')):
            continue
        # Allow common labels such as sequence=ACGT... or annotation=++--.
        if '=' in line:
            line = line.split('=', 1)[1]
        token = line.strip()
        if token and all(c.upper() in 'ACGTN' for c in token):
            seq_parts.append(token.upper())
        elif token and all(c in '+-.' for c in token):
            ann_parts.append(token.replace('.', '-'))
        else:
            # Some course files may place sequence and annotation as two tokens on one line.
            for part in raw.split():
                part = part.strip()
                if part and all(c.upper() in 'ACGTN' for c in part):
                    seq_parts.append(part.upper())
                elif part and all(c in '+-.' for c in part):
                    ann_parts.append(part.replace('.', '-'))
    seq = ''.join(seq_parts)
    ann = ''.join(ann_parts) if ann_parts else None
    if not seq:
        raise ValueError(f'No DNA sequence found in {path}; inspect the file format and adapt read_hmm_data().')
    if ann is not None and len(ann) != len(seq):
        raise ValueError(f'Sequence length ({len(seq)}) != annotation length ({len(ann)}) in {path}.')
    return seq, ann


def default_model(p_high_stay=.99, p_low_stay=.99, high_gc=.60, low_gc=.40):
    """Return transition, emission, and initial probabilities.

    The assignment text says the starter program has hard-coded values, but those values
    are not present in the attached problem statement. These editable defaults encode
    60% GC for H and 40% GC for L; replace stay probabilities with the starter's values.
    """
    trans = {'H': {'H': p_high_stay, 'L': 1-p_high_stay},
             'L': {'L': p_low_stay, 'H': 1-p_low_stay}}
    emit = {
        'H': {'A': (1-high_gc)/2, 'T': (1-high_gc)/2, 'C': high_gc/2, 'G': high_gc/2, 'N': .25},
        'L': {'A': (1-low_gc)/2, 'T': (1-low_gc)/2, 'C': low_gc/2, 'G': low_gc/2, 'N': .25},
    }
    initial = {'H': .5, 'L': .5}
    return trans, emit, initial


def viterbi(sequence, trans, emit, initial):
    """Decode most likely state path in O(length * states^2) time."""
    if not sequence:
        return ''
    neg_inf = float('-inf')
    back = [{}]
    first = sequence[0]
    for s in STATES:
        back[0][s] = None
    rows = [{s: math.log(max(initial[s], 1e-300)) + math.log(max(emit[s].get(first, emit[s].get('N', .25)), 1e-300)) for s in STATES}]
    for i in range(1, len(sequence)):
        cur, prev_back = {}, {}
        base = sequence[i]
        for s in STATES:
            options = [(rows[-1][p] + math.log(max(trans[p][s], 1e-300)), p) for p in STATES]
            best_score, best_prev = max(options)
            cur[s] = best_score + math.log(max(emit[s].get(base, emit[s].get('N', .25)), 1e-300))
            prev_back[s] = best_prev
        rows.append(cur)
        back.append(prev_back)
    state = max(STATES, key=lambda s: rows[-1][s])
    path = [state]
    for i in range(len(sequence)-1, 0, -1):
        state = back[i][state]
        path.append(state)
    path.reverse()
    return ''.join(path)


def annotation_to_states(annotation):
    return ''.join('H' if c == '+' else 'L' for c in annotation)


def states_to_annotation(states):
    return ''.join('+' if c == 'H' else '-' for c in states)


def runs(states):
    if not states:
        return []
    out, start = [], 0
    for i in range(1, len(states)):
        if states[i] != states[i-1]:
            out.append((states[start], i-start))
            start = i
    out.append((states[start], len(states)-start))
    return out


def accuracy(predicted, truth):
    if len(predicted) != len(truth):
        raise ValueError('Predicted and true annotations must have equal length.')
    return sum(a == b for a, b in zip(predicted, truth)) / len(truth) if truth else 0.0


def train_from_annotation(sequence, annotation, alpha=1.0):
    """Supervised estimate of transition/emission probabilities with Laplace smoothing."""
    states = annotation_to_states(annotation)
    trans_counts = {s: Counter({t: alpha for t in STATES}) for s in STATES}
    emit_counts = {s: Counter({b: alpha for b in BASES}) for s in STATES}
    for base, state in zip(sequence, states):
        if base in BASES:
            emit_counts[state][base] += 1
    for a, b in zip(states, states[1:]):
        trans_counts[a][b] += 1
    trans = {s: {} for s in STATES}
    for s in STATES:
        total = sum(trans_counts[s].values())
        trans[s] = {t: trans_counts[s][t] / total for t in STATES}
    emit = {s: {} for s in STATES}
    for s in STATES:
        total = sum(emit_counts[s].values())
        emit[s] = {b: emit_counts[s][b] / total for b in BASES}
        emit[s]['N'] = .25
    initial_counts = Counter({'H': alpha, 'L': alpha})
    if states:
        initial_counts[states[0]] += 1
    total = sum(initial_counts.values())
    initial = {s: initial_counts[s] / total for s in STATES}
    return trans, emit, initial


def summarize_runs(state_string):
    by_state = {'H': [], 'L': []}
    for state, length in runs(state_string):
        by_state[state].append(length)
    result = {}
    for state, lengths in by_state.items():
        result['High-GC' if state == 'H' else 'Low-GC'] = {
            'number_of_runs': len(lengths),
            'mean_duration': sum(lengths)/len(lengths) if lengths else 0,
            'duration_counts': {str(k): v for k, v in sorted(Counter(lengths).items())}
        }
    return result


def plot_durations(truth, pred, out):
    try:
        import matplotlib.pyplot as plt
    except ImportError as e:
        raise SystemExit('Plotting needs matplotlib. Install with: python -m pip install matplotlib') from e
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, state, label in zip(axes, ('H','L'), ('High-GC','Low-GC')):
        for data, name in ((truth, 'Correct'), (pred, 'Viterbi')):
            lengths = [n for s, n in runs(data) if s == state]
            if lengths:
                ax.hist(lengths, bins=30, alpha=.55, label=name)
        ax.set(title=label, xlabel='Run length (bases)', ylabel='Count')
        ax.legend()
    fig.tight_layout(); Path(out).parent.mkdir(parents=True, exist_ok=True); fig.savefig(out, dpi=160); plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('data', help='file containing DNA sequence and optional +/− annotation')
    p.add_argument('--p-high-stay', type=float, default=.99, help='P(H->H), default .99; replace with course starter value')
    p.add_argument('--p-low-stay', type=float, default=.99, help='P(L->L), default .99; replace with course starter value')
    p.add_argument('--train', help='annotated training file; estimate transition/emission probabilities from it')
    p.add_argument('--plot', help='optional duration histogram PNG')
    p.add_argument('--out', help='optional JSON output file')
    args = p.parse_args()
    seq, truth_ann = read_hmm_data(args.data)
    if args.train:
        train_seq, train_ann = read_hmm_data(args.train)
        if not train_ann:
            p.error('--train file must contain a correct +/− annotation')
        trans, emit, initial = train_from_annotation(train_seq, train_ann)
        model_source = f'supervised training file: {args.train}'
    else:
        trans, emit, initial = default_model(args.p_high_stay, args.p_low_stay)
        model_source = 'editable default parameters; check against course starter'
    predicted_states = viterbi(seq, trans, emit, initial)
    predicted_ann = states_to_annotation(predicted_states)
    result = {
        'file': args.data, 'sequence_length': len(seq), 'model_source': model_source,
        'transition_probabilities': trans,
        'expected_state_durations_if_geometric': {
            'High-GC': (1/trans['H']['L']) if trans['H']['L'] else 'infinite',
            'Low-GC': (1/trans['L']['H']) if trans['L']['H'] else 'infinite'
        },
        'predicted_accuracy': accuracy(predicted_ann, truth_ann) if truth_ann else None,
        'correct_duration_distribution': summarize_runs(annotation_to_states(truth_ann)) if truth_ann else None,
        'predicted_duration_distribution': summarize_runs(predicted_states),
        'predicted_annotation': predicted_ann
    }
    print(json.dumps(result, indent=2))
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(result, indent=2), encoding='utf-8')
    if args.plot and truth_ann:
        plot_durations(annotation_to_states(truth_ann), predicted_states, args.plot)

if __name__ == '__main__':
    main()
