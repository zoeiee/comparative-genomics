import unittest
from seqalign import needleman_wunsch, edit_distance
from dotplot import find_seed_hits, find_inversion_hits, reverse_complement
from viterbi import default_model, viterbi, states_to_annotation, summarize_runs

class TestAlgorithms(unittest.TestCase):
    def test_alignment_score(self):
        score, a, b, F, T = needleman_wunsch('GAT','GT',2,-1,-2)
        self.assertEqual(score, 2)
        self.assertEqual(len(a), len(b))

    def test_edit_distance(self):
        self.assertEqual(edit_distance('ATGC','ATGA'), 1)
        self.assertEqual(edit_distance('ABC','ABC'), 0)

    def test_exact_kmer(self):
        hits = find_seed_hits('AACCGG','CCGGTT',2,1)
        self.assertIn((2,0), hits)

    def test_reverse_complement(self):
        self.assertEqual(reverse_complement('AAGC'), 'GCTT')
        hits = find_inversion_hits('AACCGGTT','AACCGGTT',4,1)
        self.assertTrue(isinstance(hits,list))

    def test_viterbi_length(self):
        trans, emit, initial = default_model(.9,.9)
        states = viterbi('GGGGAAAACCCC',trans,emit,initial)
        self.assertEqual(len(states), 12)
        self.assertTrue(set(states) <= {'H','L'})

    def test_duration_summary(self):
        summary = summarize_runs('HHLLL')
        self.assertEqual(summary['High-GC']['mean_duration'], 2)
        self.assertEqual(summary['Low-GC']['mean_duration'], 3)

if __name__ == '__main__':
    unittest.main()
