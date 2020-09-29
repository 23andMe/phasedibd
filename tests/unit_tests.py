import numpy as np
import pandas as pd
import os
import sys
import unittest
import phasedibd as ibd

TEST_DATA_PATH = os.path.dirname(os.path.abspath(__file__)) + '/data/'


class TestPbwtBasics(unittest.TestCase):


    def test_PBWT_default_haplotypes(self):
        print("\nTesting non-templated PBWT matches over default haplotypes...")
        haplotypes = ibd.HaplotypeAlignment()
        tpbwt = ibd.TPBWTAnalysis(template=[[1]])
        ibd_segs = tpbwt.compute_ibd(haplotypes, L_m=3, L_f=0, missing_site_threshold=1, verbose=False)
        true_segs = pd.DataFrame({'chromosome': ['1'] * 3,
                                  'start': [0,0,1],
                                  'end': [3,5,4],
                                  'start_cm': [0,0,1],
                                  'end_cm': [3,5,4],
                                  'start_bp': [0,0,1],
                                  'end_bp': [3,5,4],
                                  'id1': [0,1,2],
                                  'id2': [7,6,3],
                                  'id1_haplotype': [0,0,0],
                                  'id2_haplotype': [0,0,0]})
        #print(ibd_segs)
        #print(true_segs)
        #self.assertTrue(np.all(ibd_segs.eq(true_segs)))
        self.assertTrue(np.array_equal(ibd_segs,true_segs))
    
    def test_vcf_no_map(self):
        print("\nTesting IBD compute over VCF file with out genetic map...")
        haplotypes = ibd.VcfHaplotypeAlignment(TEST_DATA_PATH + 'test.vcf')
        tpbwt = ibd.TPBWTAnalysis(template=[[1]])
        ibd_segs = tpbwt.compute_ibd(haplotypes, L_m=5, L_f=3.0, use_phase_correction=False, verbose=False)
        #print(ibd_segs)
        self.assertTrue(ibd_segs.shape[0] == 2)
   

    def test_vcf_with_map(self):
        print("\nTesting IBD compute over VCF file with genetic map...")
        haplotypes = ibd.VcfHaplotypeAlignment(TEST_DATA_PATH + 'test.vcf', 
                                               TEST_DATA_PATH + 'test.map')
        tpbwt = ibd.TPBWTAnalysis(template=[[1]])
        ibd_segs = tpbwt.compute_ibd(haplotypes, L_m=5, L_f=3.0, use_phase_correction=False, verbose=False)
        #print(ibd_segs)
        self.assertTrue(ibd_segs.shape[0] == 1)


    def test_pedigree_simulated_vcf(self):
        print("\nTesting simulated in-sample IBD compute: closely related samples in VCF simulated from 1KGP 4 generation pedigree...")
        est_ibd = pd.DataFrame()
        chromosomes = [str(x) for x in range(1, 23)]
        for chromo in chromosomes:
            vcf_path = TEST_DATA_PATH + 'pedigree_vcf/' + chromo + '.vcf'
            map_path = TEST_DATA_PATH + 'pedigree_vcf/maps/' + chromo + '.map'
            haplotypes = ibd.VcfHaplotypeAlignment(vcf_path, map_path)
            tpbwt = ibd.TPBWTAnalysis()
            ibd_segs = tpbwt.compute_ibd(haplotypes, use_phase_correction=False, verbose=False)
            est_ibd = est_ibd.append(ibd_segs)
        est_ibd['length'] = est_ibd['end_bp'] - est_ibd['start_bp']
        true_ibd = pd.read_csv(TEST_DATA_PATH + 'pedigree_vcf/true_ibd.csv')
        true_ibd = true_ibd[true_ibd['chromosome'] != 'X']
        # compute error in total IBD (in bp) and num segments per pair
        pairwise_error_ibd = []
        pairwise_error_num_seg = []
        for id1 in range(11):
            for id2 in range(id1 + 1, 12):
                m1 = (true_ibd['id1'] == id1) & (true_ibd['id2'] == id2)
                m2 = (est_ibd['id1'] == id1) & (est_ibd['id2'] == id2)
                s_true = sum(true_ibd[m1]['length'])
                s_est = sum(est_ibd[m2]['length'])
                c_true = true_ibd[m1].shape[0]
                c_est = est_ibd[m2].shape[0]
                if s_true != 0:
                    pairwise_error_ibd.append((s_est - s_true)/float(s_true))
                else:
                    pairwise_error_ibd.append(s_est - s_true)
                if c_true != 0:
                    pairwise_error_num_seg.append((c_est - c_true)/float(c_true))
                else:
                    pairwise_error_num_seg.append(c_est - c_true)
        # mean pairwise error in total IBD should be < 0.001
        self.assertTrue(abs(np.mean(pairwise_error_ibd)) < 0.001)
        # mean pairwise error in num IBD segments should be < 0.07
        self.assertTrue(abs(np.mean(pairwise_error_num_seg)) < 0.07)


    def test_1kgp_vcf(self):
        print("\nTesting 1KGP in-sample IBD compute: distantly related samples in VCF...")
        ibd_results = pd.DataFrame()
        for chromo in ['1', '2', '15', '22']:
            vcf_path = TEST_DATA_PATH + '1kgp_vcf/' + chromo + '.vcf'
            map_path = TEST_DATA_PATH + '1kgp_vcf/' + chromo + '.map'
            haplotypes = ibd.VcfHaplotypeAlignment(vcf_path, map_path)
            tpbwt = ibd.TPBWTAnalysis()
            ibd_segs = tpbwt.compute_ibd(haplotypes, L_f=4.0, use_phase_correction=False, verbose=False)
            ibd_results = ibd_results.append(ibd_segs)

        # we expect 2 IBD segments > 4 cM
        self.assertTrue(ibd_results.shape[0] == 2)
        # the segment on chromo 1 is self IBD
        id1 = ibd_results[ibd_results['chromosome'] == '1']['id1'].iloc[0]
        id2 = ibd_results[ibd_results['chromosome'] == '1']['id2'].iloc[0]
        self.assertTrue(id1 == id2)
        # the other segment should be between on chromo 15 shared between individual 3 and 7
        id1 = ibd_results[ibd_results['chromosome'] == '15']['id1'].iloc[0]
        id2 = ibd_results[ibd_results['chromosome'] == '15']['id2'].iloc[0]
        estimated_pair = set([id1, id2])
        true_pair = set([3, 7])
        self.assertTrue(estimated_pair == true_pair)


    def test_1kgp_compress(self):
        print("\nTesting 1KGP in-sample IBD compute: comparing results from TPBWT-compressed samples and VCF samples...")
        tpbwt = ibd.TPBWTAnalysis()
        uncompressed_ibd_results = pd.DataFrame()
        for chromo in ['1', '2', '15', '22']:
            vcf_path = TEST_DATA_PATH + '1kgp_vcf/' + chromo + '.vcf'
            map_path = TEST_DATA_PATH + '1kgp_vcf/' + chromo + '.map'
            haplotypes = ibd.VcfHaplotypeAlignment(vcf_path, map_path)
            tpbwt.compress_alignment('compressed_haplotypes_1kgp/', haplotypes, verbose=False)
            ibd_segs = tpbwt.compute_ibd(haplotypes, L_f=4.0, verbose=False)
            uncompressed_ibd_results = uncompressed_ibd_results.append(ibd_segs)
        haplotypes_compressed_1kgp = ibd.CompressedHaplotypeAlignment('compressed_haplotypes_1kgp/')
        compressed_ibd_results = tpbwt.compute_ibd(haplotypes_compressed_1kgp, L_f=4.0, verbose=False)
        # we expect 2 IBD segments > 4 cM
        self.assertTrue(uncompressed_ibd_results.shape[0] == 2)
        self.assertTrue(compressed_ibd_results.shape[0] == 2)
        # the segments should be identical
        self.assertTrue(sum(uncompressed_ibd_results['start_cm']) == sum(uncompressed_ibd_results['start_cm']))
        self.assertTrue(sum(uncompressed_ibd_results['end_cm']) == sum(uncompressed_ibd_results['end_cm']))


if __name__ == '__main__':
    unittest.main()
