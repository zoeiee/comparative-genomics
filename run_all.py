#!/usr/bin/env python3
"""Run available parts of MIT Computational Biology PS1 using files in ./data."""
import json
from pathlib import Path
from seqalign import read_fasta, needleman_wunsch, edit_distance
from dotplot import find_seed_hits, find_inversion_hits, plot_hits
from viterbi import read_hmm_data, default_model, viterbi, states_to_annotation, annotation_to_states, accuracy, summarize_runs, train_from_annotation, plot_durations

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data'
OUT = ROOT / 'results'
OUT.mkdir(exist_ok=True)


def first_existing(*names):
    for name in names:
        path = DATA / name
        if path.exists():
            return path
    return None


def main():
    report = {'note': 'Numbers are computed only for files present in ./data. Alignment score defaults must be set to match the provided course skeleton.'}
    # Q1: gene alignments/distances
    gene_pairs = [
        ('human_mouse_HoxA13', ('human-hoxa13.fa','human_hoxa13.fa','human-HoxA13.fa'), ('mouse-hoxa13.fa','mouse_hoxa13.fa','mouse-HoxA13.fa')),
        ('human_mouse_HoxD13', ('human-hoxd13.fa','human_hoxd13.fa','human-HoxD13.fa'), ('mouse-hoxd13.fa','mouse_hoxd13.fa','mouse-HoxD13.fa')),
        ('human_HoxA13_HoxD13', ('human-hoxa13.fa','human_hoxa13.fa','human-HoxA13.fa'), ('human-hoxd13.fa','human_hoxd13.fa','human-HoxD13.fa')),
        ('mouse_HoxA13_HoxD13', ('mouse-hoxa13.fa','mouse_hoxa13.fa','mouse-HoxA13.fa'), ('mouse-hoxd13.fa','mouse_hoxd13.fa','mouse-HoxD13.fa')),
    ]
    q1 = {}
    for label, names1, names2 in gene_pairs:
        p1, p2 = first_existing(*names1), first_existing(*names2)
        if p1 and p2:
            a, b = read_fasta(p1), read_fasta(p2)
            score, aa, bb, F, T = needleman_wunsch(a, b, match=1, mismatch=-1, gap=-1)
            dist = edit_distance(a, b)
            q1[label] = {'file1': p1.name, 'file2': p2.name, 'length1': len(a), 'length2': len(b), 'alignment_score_default_1_m1_gm1': score, 'edit_distance': dist}
            (OUT / f'{label}_alignment.txt').write_text(f'{aa}\n{bb}\nScore={score}\nEdit distance={dist}\n', encoding='utf-8')
    report['Q1_gene_comparisons'] = q1
    orth = q1.get('human_mouse_HoxA13')
    paralogs = [q1[k]['edit_distance'] for k in ('human_HoxA13_HoxD13','mouse_HoxA13_HoxD13') if k in q1]
    if orth and paralogs and orth['edit_distance'] > 0:
        avg_para = sum(paralogs)/len(paralogs)
        report['Q1_WGD_estimate_Myr'] = 70 * avg_para / orth['edit_distance']
        report['Q1_WGD_estimate_assumption'] = 'Simple proportional-distance clock; equal rates; average within-species A13/D13 paralog distance divided by human/mouse A13 ortholog distance.'
    # Q2: required seed settings
    human = first_existing('human-hoxa-region.fa','human_hoxa_region.fa')
    mouse = first_existing('mouse-hoxa-region.fa','mouse_hoxa_region.fa')
    if human and mouse:
        a, b = read_fasta(human), read_fasta(mouse)
        experiments = [('exact_30mer',30,1),('exact_100mer',100,1),('spaced_60_every2',60,2),('spaced_90_every3',90,3),('spaced_120_every4',120,4)]
        q2 = {}
        tol = max(1, int(.01*max(len(a),len(b))))
        for label, span, step in experiments:
            hits = find_seed_hits(a,b,span,step)
            near = sum(abs(x-y) <= tol for x,y in hits)
            q2[label] = {'span':span,'step':step,'sampled_bases':len(range(0,span,step)),'hits':len(hits),'near_diagonal_tolerance':tol,'near_diagonal_hits':near,'near_diagonal_percent':100*near/len(hits) if hits else 0}
            plot_hits(hits,len(a),len(b),OUT/f'{label}.png',label)
        mod = first_existing('human-hoxa-region-modified.fa','human_hoxa_region_modified.fa')
        if mod:
            invseq = read_fasta(mod)
            forward = find_seed_hits(invseq,b,30,1)
            reverse = find_inversion_hits(invseq,b,30,1)
            plot_hits(forward,len(invseq),len(b),OUT/'inversion_dotplot.png','Artificial inversion: 30-mer seeds',reverse)
            q2['inversion_demo'] = {'forward_hits':len(forward),'reverse_complement_hits':len(reverse),'plot':'results/inversion_dotplot.png'}
        report['Q2_dotplots'] = q2
    # Q3: HMM data; train on hmmgen if available, evaluate all present files
    hmm_files = [p for p in [first_existing('hmmgen'), first_existing('mystery1'), first_existing('mystery2'), first_existing('mystery3')] if p]
    if hmm_files:
        train_pair = None
        train_path = first_existing('hmmgen')
        if train_path:
            train_pair = read_hmm_data(train_path)
        q3 = {}
        for path in hmm_files:
            seq, ann = read_hmm_data(path)
            if train_pair and train_path != path and train_pair[1]:
                trans, emit, initial = train_from_annotation(*train_pair)
                source = 'trained on hmmgen annotations'
            else:
                trans, emit, initial = default_model(.99,.99)
                source = 'default .99 stay probabilities; replace with original starter values'
            pred_states = viterbi(seq,trans,emit,initial)
            pred_ann = states_to_annotation(pred_states)
            q3[path.name] = {'length':len(seq),'model_source':source,'accuracy':accuracy(pred_ann,ann) if ann else None,'correct_durations':summarize_runs(annotation_to_states(ann)) if ann else None,'predicted_durations':summarize_runs(pred_states)}
            if ann:
                plot_durations(annotation_to_states(ann),pred_states,OUT/f'{path.name}_durations.png')
        report['Q3_HMM'] = q3
    (OUT/'summary.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))
    print(f'\nSaved report and outputs in: {OUT}')
    print('Missing input files are skipped. See README.md for required filenames and model caveats.')

if __name__ == '__main__':
    main()
