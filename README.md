# MIT Computational Biology PS1 — runnable Python project

This project implements the main computational work in the attached MIT 6.047/6.878/HST.507 Problem Set 1: global sequence alignment and evolutionary distances (Q1), k-mer/spaced-seed dotplots and inversion matching (Q2), and a two-state GC-content HMM with Viterbi decoding (Q3). Q4 is a personal student profile, so a truthful fill-in template is provided instead of invented personal history.

## Important data/scoring caveat

The uploaded file contains the problem statement, not the original starter scripts or sequence/HMM data. Put the official files in `data/` (see `data/README.md`). The problem statement says `ps1-seqalign.py` and `viterbi.py` have specified scoring/model parameters, but those parameter values are not visible in the uploaded text. This implementation exposes the alignment scores and HMM self-transition probabilities as command-line options; for exact assignment numbers, set them to the values in the official starter files. Do not report example/default outputs as the official answer.

## Setup (macOS/Linux)

```bash
cd mit_compbio_ps1
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## Q1 — Needleman–Wunsch and distances

```bash
python seqalign.py data/seq1.fa data/seq2.fa --match 1 --mismatch -1 --gap -1 --matrix-out results/alignment_matrix.csv --distance
```

- `--match`, `--mismatch`, and `--gap` must be set to the values in the original course skeleton for the required score.
- `--matrix-out` writes both the score and traceback matrices.
- `--distance` calculates Levenshtein edit distance, a nonnegative dissimilarity where identical sequences have distance zero.
- The short example can be run by creating two FASTA files containing `CTAAGTACT` and `CATTA` and passing the matching starter scoring parameters.

For a quick comparison of HoxA13/HoxD13 pairs, place the four gene FASTA files in `data/` and run `python run_all.py`. It reports ortholog and within-species paralog edit distances and a simple proportional-clock estimate. This estimate assumes equal rates and is not a model-corrected evolutionary distance.

## Q2 — dotplots and hashing

```bash
# Exact matching 30-mers
python dotplot.py data/human-hoxa-region.fa data/mouse-hoxa-region.fa results/30mer.png --span 30 --step 1 --summary-out results/30mer.json

# Exact 100-mers
python dotplot.py data/human-hoxa-region.fa data/mouse-hoxa-region.fa results/100mer.png --span 100 --step 1

# 60-mer matching every other base (30 sampled bases)
python dotplot.py data/human-hoxa-region.fa data/mouse-hoxa-region.fa results/spaced60.png --span 60 --step 2

# 90-mer matching every third base
python dotplot.py data/human-hoxa-region.fa data/mouse-hoxa-region.fa results/spaced90.png --span 90 --step 3

# 120-mer matching every fourth base
python dotplot.py data/human-hoxa-region.fa data/mouse-hoxa-region.fa results/spaced120.png --span 120 --step 4

# Artificial inversion; displays forward and reverse-complement hits
python dotplot.py data/human-hoxa-region-modified.fa data/mouse-hoxa-region.fa results/inversion.png --span 30 --step 1 --inversions

# Q2(b)(v): 100-mer pairs with at most two mismatches in every 6-base block
python dotplot.py data/human-hoxa-region.fa data/mouse-hoxa-region.fa results/local-rule.png --local-mismatch-rule
```

`step=2/3/4` samples every second/third/fourth base starting at offset 0. This is one consistent interpretation of the assignment's spaced-seed wording. The near-diagonal threshold defaults to 1% of the longer sequence length and is printed in the JSON; if the official helper uses a different threshold, change `--diag-tolerance` so your numbers are comparable. Q2(b)(v) uses seed-and-verify candidate generation and can be computationally expensive on highly repetitive data.

## Q3 — HMM/Viterbi

```bash
# Decode with default model; change stay probabilities to match official viterbi.py
python viterbi.py data/hmmgen --p-high-stay 0.99 --p-low-stay 0.99 --out results/hmmgen.json --plot results/hmmgen_durations.png

# Fit transitions/emissions from annotated hmmgen, then decode a mystery sequence
python viterbi.py data/mystery1 --train data/hmmgen --out results/mystery1.json --plot results/mystery1_durations.png
python viterbi.py data/mystery2 --train data/hmmgen --out results/mystery2.json
python viterbi.py data/mystery3 --train data/hmmgen --out results/mystery3.json
```

The parser expects DNA letters (`A/C/G/T/N`) and an annotation made of `+` (High-GC) and `-` (Low-GC); they can be on separate lines. If your original data file has another format, adapt `read_hmm_data()` in `viterbi.py`. The provided problem statement does not include the actual hard-coded transition/emission values, so the defaults are illustrative until checked against the original starter. Training on `hmmgen` is supervised and avoids using the mystery sequence's own annotation to train its model.

For a state with self-transition probability `a_kk`, duration is geometric:

`P(D=d) = (a_kk ** (d-1)) * (1-a_kk)`, for `d=1,2,...`

`E[D] = 1/(1-a_kk)`.

## Run all files present

```bash
python run_all.py
```

This checks for expected filenames, runs available comparisons/plots, and writes `results/summary.json`. Missing data files are skipped and reported in `data/README.md`; numerical results cannot be produced honestly until the official inputs are present.

## Tests

```bash
python -m unittest discover -s tests -v
```

## What this project covers / what still needs human input

- Q1(a–d): code, alignments, score matrices, edit distances, and a simple duplication-date estimate.
- Q2(a–e): exact/spaced seeds, hit counts, near-diagonal fraction, plots, local mismatch rule, and reverse-complement/inversion visualization.
- Q3(a–e): geometric duration formula, Viterbi decoding, accuracy, run-length distributions, supervised parameter estimation, and plots. The written explanation of GENSCAN's non-geometric duration modeling belongs in the report; see your course paper/lecture notes. The code cannot infer the course's missing starter parameters.
- Q4: `student_profile_template.md`; fill in your own actual skills, experience, interests, and project preferences.
