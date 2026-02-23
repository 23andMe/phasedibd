# cython: profile=True
# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
import itertools
import math
import numpy as np
import os
import pandas as pd
import sys
import timeit
cimport numpy as np
from libc.stdio cimport *
from libc.stdlib cimport malloc, free, realloc
from libc.string cimport *
from libcpp cimport bool
from posix.unistd cimport access, F_OK


cdef class TPBWT:
    """
    Class that contains the core functionality to:
    - compute IBD segments among haplotypes using the TPBWT algorithm
    - write TPBWT-compressed haplotypes to binary files
    - write compressed haplotype metadata to files

    Setting template to [[1]] collapses the TPBWT down to the PBWT.
    """

    def __init__(self, template=None):
        
        self.num_all_x_haplotypes = 0
        self.num_all_y_haplotypes = 0

        # set up the template function
        if template is None:
            self.template = [[1, 0, 1, 0],
                             [0, 1, 0, 1],
                             [1, 1, 0, 0],
                             [0, 0, 1, 1],
                             [1, 0, 0, 1],
                             [0, 1, 1, 0]]
        else:
            if type(template) is not list or \
                    len(template) == 0 or \
                    type(template[0]) is not list or \
                    len(template[0]) == 0:
                raise Exception('Invalid template.')
            self.template = template
        self.num_templates = len(self.template)
        self.template_length = len(self.template[0])
        self.template_ptr = <uint8_t **>malloc(self.num_templates * sizeof(uint8_t*))
        for i in xrange(self.num_templates):
            self.template_ptr[i] = <uint8_t *>malloc(self.template_length * sizeof(uint8_t))
            for j in xrange(self.template_length):
                if self.template[i][j] == 0:
                    self.template_ptr[i][j] = 0
                else:
                    self.template_ptr[i][j] = 1

        # allocate empty arrays to hold matching segments
        self.num_segments = 0
        self.segments_capacity = 10000000
        self.segments_id1 = <long *>malloc(self.segments_capacity * sizeof(long))
        self.segments_id2 = <long *>malloc(self.segments_capacity * sizeof(long))
        self.segments_id1_haplotype = <uint8_t *>malloc(self.segments_capacity * sizeof(uint8_t))
        self.segments_id2_haplotype = <uint8_t *>malloc(self.segments_capacity * sizeof(uint8_t))
        self.segments_start = <uint32_t *>malloc(self.segments_capacity * sizeof(uint32_t))
        self.segments_end = <uint32_t *>malloc(self.segments_capacity * sizeof(uint32_t))
        self.segments_start_cm = <double *>malloc(self.segments_capacity * sizeof(double))
        self.segments_end_cm = <double *>malloc(self.segments_capacity * sizeof(double))
        self.segments_start_bp = <uint32_t *>malloc(self.segments_capacity * sizeof(uint32_t))
        self.segments_end_bp = <uint32_t *>malloc(self.segments_capacity * sizeof(uint32_t))
        self.segments_chromosome = []
        self.segment_output_file = NULL

        self.current_matches_start = NULL
        self.current_matches_end = NULL
        self.current_switch_error = NULL


    def __dealloc__(self):

        free(self.segments_id1)
        free(self.segments_id2)
        free(self.segments_id1_haplotype)
        free(self.segments_id2_haplotype)
        free(self.segments_start)
        free(self.segments_end)
        free(self.segments_start_cm)
        free(self.segments_end_cm)
        free(self.segments_start_bp)
        free(self.segments_end_bp)
        cdef uint32_t num_cols
        if self.current_matches_start != NULL:
            num_cols = <uint32_t>(sizeof(self.current_matches_start)/sizeof(self.current_matches_start[0]))
            for i in xrange(num_cols):
                free(self.current_matches_start[i])
                free(self.current_matches_end[i])
            free(self.current_matches_start)
            free(self.current_matches_end)
            free(self.current_switch_error)


    cdef uint8_t write_haplotype_data(self, path):
        """
        Writes haplotype alignment metadata for .tpbwt files.
        These files are used when reading the compressed haplotypes.
        Writes the number of all haplotypes, and then the haplotype names and numbers.
        Then write the number of sites, the chromosome, and the genetic and
        physical positions of each sites.
        """
        # total num haplotypes
        cdef uint32_t num_all_x_haplotypes = self.num_all_x_haplotypes
        cdef uint32_t num_all_y_haplotypes = self.num_all_y_haplotypes
        cdef uint32_t num_all_haplotypes = num_all_x_haplotypes + num_all_y_haplotypes
        
        # write to file
        cdef uint32_t i 
        cdef long temp
        cdef uint8_t temp2
        cdef FILE *f = fopen(path, "wb") 
        fwrite(&num_all_haplotypes, sizeof(uint32_t), 1, f)
        for i in xrange(num_all_x_haplotypes):
            temp = self.x_haplotypes.get_haplotype_name(i)
            fwrite(&temp, sizeof(long), 1, f)
            temp2 = self.x_haplotypes.get_haplotype_number(i)
            fwrite(&temp2, sizeof(uint8_t), 1, f)
        for i in xrange(num_all_y_haplotypes):
            temp = self.y_haplotypes.get_haplotype_name(i)
            fwrite(&temp, sizeof(long), 1, f)
            temp2 = self.y_haplotypes.get_haplotype_number(i)
            fwrite(&temp2, sizeof(uint8_t), 1, f)
        cdef uint32_t num_sites = self.x_haplotypes.get_num_sites()
        fwrite(&num_sites, sizeof(uint32_t), 1, f)
        cdef uint8_t chromo_length = len(self.chromosome)
        fwrite(&chromo_length, sizeof(uint8_t), 1, f)
        chromo_temp = self.chromosome.encode()
        fwrite(<char*>chromo_temp, sizeof(char), chromo_length, f)
        cdef double gen_pos
        cdef uint32_t phy_pos
        for i in xrange(num_sites):
            gen_pos = self.x_haplotypes.get_genetic_position_from_site(i)
            fwrite(&gen_pos, sizeof(double), 1, f)
            phy_pos = self.x_haplotypes.get_physical_position_from_site(i)
            fwrite(&phy_pos, sizeof(uint32_t), 1, f)
        fclose(f)
    
    
    
    cdef uint8_t allocate_current_match_arrays(self):
        """
        In the case of an in-sample analysis num_all_y_haplotypes = 0 and the start/end
        arrays should be num_all_x_haplotypes by num_all_x_haplotypes in size.
        In the case of an out-of-sample analysis num_all_y_haplotypes > 0 and the start/end
        arrays should be num_all_x_haplotypes by num_all_y_haplotypes in size.

        The switch error array should always be num_all_x_haplotypes + num_all_y_haplotypes
        in length.
        """
        cdef uint32_t i
        cdef uint32_t j
        cdef uint32_t num_rows
        cdef uint32_t num_cols

        # first free any already allocated memory
        if self.current_matches_start != NULL:
            num_rows = <uint32_t>(sizeof(self.current_matches_start)/sizeof(self.current_matches_start[0]))
            for i in xrange(num_rows):
                free(self.current_matches_start[i])
                free(self.current_matches_end[i])
            free(self.current_matches_start)
            free(self.current_matches_end)
            free(self.current_switch_error)

        # now initialize the arrays with the correct size
        num_rows = self.num_all_x_haplotypes
        num_cols = self.num_all_y_haplotypes
        if self.num_all_y_haplotypes == 0:
            num_cols = self.num_all_x_haplotypes
        self.current_matches_start = <uint32_t **>malloc(num_rows * sizeof(uint32_t*))
        self.current_matches_end = <uint32_t **>malloc(num_rows * sizeof(uint32_t*))
        
        for i in xrange(num_rows):
            self.current_matches_start[i] = <uint32_t *>malloc(num_cols * sizeof(uint32_t))
            self.current_matches_end[i] = <uint32_t *>malloc(num_cols * sizeof(uint32_t))
            for j in xrange(num_cols):
                self.current_matches_start[i][j] = 0
                self.current_matches_end[i][j] = 0
        
        cdef uint32_t num_all_haplotypes = self.num_all_x_haplotypes + self.num_all_y_haplotypes
        self.current_switch_error = <uint8_t *>malloc(num_all_haplotypes * sizeof(uint8_t))
        for i in xrange(num_all_haplotypes):
            self.current_switch_error[i] = 0


    cdef transform(self, HaplotypeAlignment x_haplotypes, HaplotypeAlignment y_haplotypes=None, char* compressed_out_path=NULL, uint32_t L_m=300, double L_f=3.0, uint32_t missing_site_threshold=10, bint use_phase_correction=True, char* segments_out_path=NULL, bint verbose=True, bint compute_ibd=True):
        """
        Method that implements the core TPBWT algorithm.
        Optionally computes IBD segments among haplotypes and optionally writes TPBWT-compressed haplotypes to binary files.
        """
        # file for output segments
        if segments_out_path != NULL:
            self.segment_output_file = fopen(segments_out_path, 'w') 
            row = 'id1,id2,id1_hap,id2_hap,start,end,start_cm,end_cm,start_bp,end_bp,chromosome\n'
            row_len = len(row)
            temp_row = row.encode()
            fwrite(<char*>temp_row, sizeof(char), row_len, self.segment_output_file)

        # setup access to the haplotype alignment
        self.x_haplotypes = x_haplotypes
        self.y_haplotypes = y_haplotypes
        self.chromosome = x_haplotypes.get_current_chromosome()
       
        cdef uint32_t num_sites_x = x_haplotypes.get_num_sites()
        cdef uint32_t num_sites_y = 0
        if y_haplotypes is not None:
            num_sites_y = y_haplotypes.get_num_sites()
            if num_sites_y != num_sites_x:
                raise Exception('The number of sites in the two haplotype alignments do not match!')
        
        cdef uint32_t num_active_x_haplotypes = x_haplotypes.get_num_active_haplotypes()
        cdef uint32_t num_all_x_haplotypes = x_haplotypes.get_num_all_haplotypes()
        cdef uint32_t num_active_y_haplotypes = 0
        cdef uint32_t num_all_y_haplotypes = 0
        if y_haplotypes is not None:
            num_active_y_haplotypes = y_haplotypes.get_num_active_haplotypes()
            num_all_y_haplotypes = y_haplotypes.get_num_all_haplotypes()
        
        # num haplotypes to be processed in this TPBWT out of the total number of haplotypes
        cdef uint32_t total_num_active_haplotypes = num_active_x_haplotypes + num_active_y_haplotypes
       
        # calculate the total number of haplotypes in both alignments
        cdef uint32_t total_num_all_haplotypes = num_all_x_haplotypes + num_all_y_haplotypes
        
        cdef Py_ssize_t i 
        cdef Py_ssize_t j

        if compute_ibd and verbose:
            if self.y_haplotypes is None:
                print('Running in-sample TPBWT on ' + str(num_sites_x) + ' sites.')
                print('Number of haplotypes: ' + str(num_all_x_haplotypes))
                print('Minimum memory required: ' + str((int(num_all_x_haplotypes) \
                                                         * int(num_all_x_haplotypes) \
                                                         * 4 * 2)/1e+9) + ' Gb')
                sys.stdout.flush()
            else:
                print('Running out-of-sample TPBWT on ' + str(num_sites_x) + ' sites.')
                print('Number of haplotypes in x: ' + str(num_all_x_haplotypes))
                print('Number of haplotypes in y: ' + str(num_all_y_haplotypes))
                print('Minimum memory required: ' + str((int(num_all_x_haplotypes) \
                                                         * int(num_all_y_haplotypes) \
                                                         * 4 * 2)/1e+9) + ' Gb')
                sys.stdout.flush()
        
        # allocate/resize current match count array if necessary
        if self.num_all_x_haplotypes != num_all_x_haplotypes \
           or self.num_all_y_haplotypes != num_all_y_haplotypes:
            self.num_all_x_haplotypes = num_all_x_haplotypes
            self.num_all_y_haplotypes = num_all_y_haplotypes
            if compute_ibd:
                self.allocate_current_match_arrays()
        elif compute_ibd and self.current_matches_start == NULL:
            self.allocate_current_match_arrays()
            
        cdef bint in_sample_analysis = num_all_y_haplotypes == 0
        cdef bint match_long_enough
        cdef bint out_of_sample_match
      
        # initialise positional prefix array
        cdef long* ppa_x = x_haplotypes.get_active_haplotype_indices()
        cdef long* ppa_y = NULL
        if y_haplotypes is not None:
            ppa_y = y_haplotypes.get_active_haplotype_indices()
        cdef long** ppa = <long **>malloc(self.num_templates * sizeof(long*))
        for i in xrange(self.num_templates):
            ppa[i] = <long *>malloc(total_num_active_haplotypes * sizeof(long))
            for j in xrange(total_num_active_haplotypes):
                if j < num_active_x_haplotypes:
                    ppa[i][j] = ppa_x[j]
                else:
                    ppa[i][j] = ppa_y[j - num_active_x_haplotypes] + num_all_x_haplotypes
        cdef long* ppa_0 = <long *>malloc(total_num_active_haplotypes * sizeof(long))
        cdef long* ppa_1 = <long *>malloc(total_num_active_haplotypes * sizeof(long))
        cdef Py_ssize_t ppa_0_i = 0
        cdef Py_ssize_t ppa_1_i = 0
       
        # write header info for TPBWT compressed haplotypes
        cdef FILE *f_out 
        if compressed_out_path != NULL:
            f_out = fopen(compressed_out_path, 'wb') 
            # write the number of haplotypes in the TPBWT, the num of active haplotype indices used,
            # the total number of haplotypes, and the number of sites in this TPBWT
            fwrite(&total_num_active_haplotypes, sizeof(uint32_t), 1, f_out)
            fwrite(&ppa[0][0], sizeof(long), total_num_active_haplotypes, f_out)
            fwrite(&total_num_all_haplotypes, sizeof(uint32_t), 1, f_out)
            fwrite(&num_sites_x, sizeof(uint32_t), 1, f_out)
        
        # initialize array to hold counts of adjacent missing sites
        cdef uint8_t* missing = <uint8_t *>malloc(total_num_all_haplotypes * sizeof(uint8_t))
        for i in xrange(total_num_all_haplotypes):
            missing[i] = 0
        
        # initialise divergence array
        cdef uint32_t** div = <uint32_t **>malloc(self.num_templates * sizeof(uint32_t*))
        for i in xrange(self.num_templates):
            div[i] = <uint32_t *>malloc(total_num_active_haplotypes * sizeof(uint32_t))
            for j in xrange(total_num_active_haplotypes):
                div[i][j] = 0
        cdef uint32_t* div_0 = <uint32_t *>malloc(total_num_active_haplotypes * sizeof(uint32_t))
        cdef uint32_t* div_1 = <uint32_t *>malloc(total_num_active_haplotypes * sizeof(uint32_t))
        cdef uint32_t div_0_i = 0
        cdef uint32_t div_1_i = 0
        
        # initialise array to hold the current allele
        cdef long* alleles = <long *>malloc(total_num_active_haplotypes * sizeof(long))

        # counts of the number of matches
        cdef Py_ssize_t num_matches_end_0 = 0
        cdef Py_ssize_t num_matches_end_1 = 0

        cdef uint32_t start = 0
        cdef uint32_t start_0 = 0
        cdef uint32_t start_1 = 0
        cdef uint32_t end = 0
        
        cdef uint32_t current_site
        cdef uint32_t last_site
        cdef long haplotype_index
        cdef long index_0
        cdef long index_1
        cdef Py_ssize_t l
        cdef Py_ssize_t u
        cdef Py_ssize_t y
        cdef uint8_t allele
        cdef uint8_t last_allele
        cdef uint8_t template_pos = 0
        cdef uint8_t t = 0
        
        cdef uint16_t current_allele_count 
       
        # Here we do some ugly setup to inline the logic to get alleles from the haplotype 
        # alignment objects. If we do not inline this logic, each time
        # x_haplotypes.get_allele(haplo, site) is called the pointer to the cdef class method will 
        # be looked up in the __pyx_vtab table which is slow. So if possible we provide explicit
        # logic here for each different type of HaplotypeAlignment.
        # Note this is because the cython compiler will not inline class methods (though it will
        # inline module level functions).
        cdef uint8_t** x_alignment
        cdef uint8_t** y_alignment
        cdef FILE *x_compressed_tpbwt_file
        cdef uint8_t x_current_allele
        cdef uint16_t x_current_allele_total
        cdef uint16_t x_current_allele_count
        cdef uint32_t x_haplotype_count
        cdef FILE *y_compressed_tpbwt_file
        cdef uint8_t y_current_allele
        cdef uint16_t y_current_allele_total
        cdef uint16_t y_current_allele_count
        cdef uint32_t y_haplotype_count
        cdef FILE *vcf_file_x
        cdef FILE *vcf_file_y
        cdef size_t vcf_sample_start_x
        cdef size_t vcf_sample_start_y
        cdef size_t vcf_num_read_lines = 0
        cdef size_t vcf_current_sample = 0
        cdef char *vcf_line = NULL
        cdef char *vcf_word
        vcf_delimiter_py = '\t'.encode()
        cdef char *vcf_delimiter = vcf_delimiter_py
        cdef size_t vcf_n = 0
        cdef ssize_t vcf_line_len = 0
        cdef size_t vcf_word_num = 0

        if type(x_haplotypes) is CompressedHaplotypeAlignment:
            x_compressed_tpbwt_file = (<CompressedHaplotypeAlignment>x_haplotypes).get_tpbwt_file()
            x_haplotype_count = 0
        elif type(x_haplotypes) is VcfHaplotypeAlignment:
            x_alignment = x_haplotypes.get_alignment()
            # start reading vcf
            vcf_file_x = fopen(x_haplotypes.get_vcf_file(), 'r')
            vcf_num_read_lines = 0
            while vcf_num_read_lines <= x_haplotypes.get_header_length():
                vcf_line_len = getline(&vcf_line, &vcf_n, vcf_file_x)
                vcf_num_read_lines += 1
            vcf_sample_start_x = x_haplotypes.get_sample_start()
        if type(y_haplotypes) is CompressedHaplotypeAlignment:
            y_compressed_tpbwt_file = (<CompressedHaplotypeAlignment>y_haplotypes).get_tpbwt_file()
            y_haplotype_count = 0
        elif type(y_haplotypes) is VcfHaplotypeAlignment:
            y_alignment = y_haplotypes.get_alignment()
            # start reading vcf
            vcf_file_y = fopen(y_haplotypes.get_vcf_file(), 'r')
            vcf_num_read_lines = 0
            while vcf_num_read_lines <= y_haplotypes.get_header_length():
                vcf_line_len = getline(&vcf_line, &vcf_n, vcf_file_y)
                vcf_num_read_lines += 1
            vcf_sample_start_y = y_haplotypes.get_sample_start()

        # iterate over variants
        for current_site in xrange(num_sites_x):
            
            # read another line of the vcf
            if type(x_haplotypes) is VcfHaplotypeAlignment:
                #x_haplotypes.update_position(current_site)
                vcf_line_len = getline(&vcf_line, &vcf_n, vcf_file_x)
                vcf_word_num = 0
                vcf_current_sample = 0
                if vcf_line_len != -1: 
                    # parse line
                    vcf_word = strtok(vcf_line, vcf_delimiter)
                    while vcf_word != NULL:
                        if vcf_word_num >= vcf_sample_start_x:
                            # get haplotype alleles for this site
                            # 48 is ascii for '0'
                            # 49 is ascii for '1'
                            if vcf_word[0] == 48:
                                x_alignment[vcf_current_sample * 2][0] = 0
                            elif vcf_word[0] == 49:
                                x_alignment[vcf_current_sample * 2][0] = 1
                            else:
                                x_alignment[vcf_current_sample * 2][0] = 2
                            if vcf_word[2] == 48:
                                x_alignment[(vcf_current_sample * 2) + 1][0] = 0
                            elif vcf_word[2] == 49:
                                x_alignment[(vcf_current_sample * 2) + 1][0] = 1
                            else:
                                x_alignment[(vcf_current_sample * 2) + 1][0] = 2
                            vcf_current_sample += 1
                        vcf_word = strtok(NULL, vcf_delimiter)
                        vcf_word_num += 1
                else:
                    fclose(vcf_file_x)
            if type(y_haplotypes) is VcfHaplotypeAlignment:
                #y_haplotypes.update_position(current_site)
                vcf_line_len = getline(&vcf_line, &vcf_n, vcf_file_y)
                vcf_word_num = 0
                vcf_current_sample = 0
                if vcf_line_len != -1: 
                    # parse line
                    vcf_word = strtok(vcf_line, vcf_delimiter)
                    while vcf_word != NULL:
                        if vcf_word_num >= vcf_sample_start_y:
                            # get haplotype alleles for this site
                            # 48 is ascii for '0'
                            # 49 is ascii for '1'
                            if vcf_word[0] == 48:
                                y_alignment[vcf_current_sample * 2][0] = 0
                            elif vcf_word[0] == 49:
                                y_alignment[vcf_current_sample * 2][0] = 1
                            else:
                                y_alignment[vcf_current_sample * 2][0] = 2
                            if vcf_word[2] == 48:
                                y_alignment[(vcf_current_sample * 2) + 1][0] = 0
                            elif vcf_word[2] == 49:
                                y_alignment[(vcf_current_sample * 2) + 1][0] = 1
                            else:
                                y_alignment[(vcf_current_sample * 2) + 1][0] = 2
                            vcf_current_sample += 1
                        vcf_word = strtok(NULL, vcf_delimiter)
                        vcf_word_num += 1
                else:
                    fclose(vcf_file_y)
            
            # iterate over templates
            for t in xrange(self.num_templates):

                # determine which templates to process for this site:
                if self.template_ptr[t][template_pos] == 1:
                    
                    # setup intermediates
                    ppa_0_i = 0
                    ppa_1_i = 0
                    div_0_i = 0
                    div_1_i = 0
                    num_matches_end_0 = 0
                    num_matches_end_1 = 0
                    current_allele_count = 0
                   
                    # find the next templated site
                    i = template_pos
                    start_0 = start_1 = current_site
                    while True:
                        i += 1
                        start_0 += 1
                        start_1 += 1
                        if i == self.template_length:
                            i = 0
                        if self.template_ptr[t][i] == 1:
                            break
                    if start_0 >= num_sites_x:
                        start_0 = start_1 = num_sites_x - 1
                    
                    # find the previous templated site
                    i = template_pos
                    last_site = current_site
                    while last_site > 0:
                        if i == 0:
                            i = self.template_length
                        i -= 1
                        last_site -= 1
                        if self.template_ptr[t][i] == 1:
                            break
                    
                    # iterate over haplotypes in reverse prefix sorted order
                    for i in xrange(total_num_active_haplotypes):

                        haplotype_index = ppa[t][i]
                        
                        # reached break in reportable matches so output current set of valid matches
                        end = last_site
                        if compute_ibd and end > div[t][i] and end - div[t][i] < L_m:
                            for u in xrange(i - (num_matches_end_0 + num_matches_end_1), i):
                                start = 0
                                for y in xrange(u + 1, i):
                                    if div[t][y] > start:
                                        start = div[t][y]
                                    match_long_enough = end > start and end - start >= L_m
                                    if match_long_enough:
                                        if alleles[u] != alleles[y]:
                                            out_of_sample_match = False
                                            if not in_sample_analysis:
                                                out_of_sample_match = (((ppa[t][u] < num_all_x_haplotypes) and \
                                                                        (ppa[t][y] >= num_all_x_haplotypes)) or \
                                                                       ((ppa[t][y] < num_all_x_haplotypes) and \
                                                                        (ppa[t][u] >= num_all_x_haplotypes)))
                                            if in_sample_analysis or out_of_sample_match:
                                                self.report_templated_segment(ppa[t][u], ppa[t][y], 
                                                        start, end, L_m, L_f, use_phase_correction)
                                    else:
                                        break
                            num_matches_end_0 = 0
                            num_matches_end_1 = 0
                        
                        # update div intermediates
                        if div[t][i] > start_0:
                            start_0 = div[t][i]
                        if div[t][i] > start_1:
                            start_1 = div[t][i]
                        
                        # get allele for current haplotype at current site
                        # for speed we need ugly inline logic here for each different type of 
                        # HaplotypeAlignment (to avoid __pyx_vtab look up)
                        # Note this is because the cython compiler will not inline class 
                        # methods (though it will inline module level functions).
                        if haplotype_index < num_all_x_haplotypes:
                            if type(x_haplotypes) is VcfHaplotypeAlignment:
                                allele = x_alignment[haplotype_index][0]
                                if allele > 1:
                                    allele = 2
                            elif type(x_haplotypes) is CompressedHaplotypeAlignment:
                                if x_haplotype_count == 0:
                                    # start reading in haplotypes for a new column of the alignment
                                    fread(&x_current_allele, sizeof(uint8_t), 1, 
                                            x_compressed_tpbwt_file)
                                    fread(&x_current_allele_total, sizeof(uint16_t), 1, 
                                            x_compressed_tpbwt_file)
                                    # keep track of the number of haplotypes in this allele run
                                    x_current_allele_count = 1
                                    # keep track of the number of haplotypes seen in this column
                                    x_haplotype_count = 1
                                elif x_current_allele_count == x_current_allele_total:
                                    # read in the next allele run for this column
                                    fread(&x_current_allele, sizeof(uint8_t), 1, 
                                            x_compressed_tpbwt_file)
                                    fread(&x_current_allele_total, sizeof(uint16_t), 1, 
                                            x_compressed_tpbwt_file)
                                    x_current_allele_count = 1
                                    x_haplotype_count += 1
                                else:
                                    # continue in the current allele run
                                    x_current_allele_count += 1
                                    x_haplotype_count += 1
                                if x_haplotype_count == num_active_x_haplotypes:
                                    # we've reached the end of this column of the alignment
                                    x_haplotype_count = 0
                                allele = x_current_allele
                            else:
                                allele = x_haplotypes.get_allele(haplotype_index, current_site)
                        else:
                            if type(y_haplotypes) is VcfHaplotypeAlignment:
                                allele = y_alignment[haplotype_index - num_all_x_haplotypes][0]
                                if allele > 1:
                                    allele = 2
                            elif type(y_haplotypes) is CompressedHaplotypeAlignment:
                                if y_haplotype_count == 0:
                                    # start reading in haplotypes for a new column of the alignment
                                    fread(&y_current_allele, sizeof(uint8_t), 1, 
                                            y_compressed_tpbwt_file)
                                    fread(&y_current_allele_total, sizeof(uint16_t), 1, 
                                            y_compressed_tpbwt_file)
                                    # keep track of the number of haplotypes in this allele run
                                    y_current_allele_count = 1
                                    # keep track of the number of haplotypes seen in this column
                                    y_haplotype_count = 1
                                elif y_current_allele_count == y_current_allele_total:
                                    # read in the next allele run for this column
                                    fread(&y_current_allele, sizeof(uint8_t), 1, 
                                            y_compressed_tpbwt_file)
                                    fread(&y_current_allele_total, sizeof(uint16_t), 1, 
                                            y_compressed_tpbwt_file)
                                    y_current_allele_count = 1
                                    y_haplotype_count += 1
                                else:
                                    # continue in the current allele run
                                    y_current_allele_count += 1
                                    y_haplotype_count += 1
                                if y_haplotype_count == num_active_y_haplotypes:
                                    # we've reached the end of this column of the alignment
                                    y_haplotype_count = 0
                                allele = y_current_allele
                            else:
                                allele = y_haplotypes.get_allele(haplotype_index - num_all_x_haplotypes, current_site)
                        # write compressed haplotypes to file
                        # note: missing/no calls are represented as 2
                        if compressed_out_path != NULL:
                            if i == total_num_active_haplotypes - 1:
                                # if we have finished with this site
                                if allele != last_allele:
                                    fwrite(&last_allele, sizeof(uint8_t), 1, f_out)
                                    fwrite(&current_allele_count, sizeof(uint16_t), 1, f_out)
                                    current_allele_count = 1
                                else:
                                    current_allele_count += 1
                                fwrite(&allele, sizeof(uint8_t), 1, f_out)
                                fwrite(&current_allele_count, sizeof(uint16_t), 1, f_out)
                            elif allele != last_allele and i != 0:
                                # write the last allele and its count
                                fwrite(&last_allele, sizeof(uint8_t), 1, f_out)
                                fwrite(&current_allele_count, sizeof(uint16_t), 1, f_out)
                                current_allele_count = 1
                            else:
                                if current_allele_count == 32767:
                                    # avoid integer overflow
                                    fwrite(&allele, sizeof(uint8_t), 1, f_out)
                                    fwrite(&current_allele_count, sizeof(uint16_t), 1, f_out)
                                    current_allele_count = 1
                                else:
                                    current_allele_count += 1

                        # handle missing/no calls
                        if allele == 2:
                            if missing[haplotype_index] < missing_site_threshold:
                                # we have not exceeded the threshold number
                                # of adjacent missing sites so extend the match
                                allele = last_allele
                            else:
                                # exceeded threshold so terminate the match
                                if last_allele == 0:
                                    allele = 1
                                else:
                                    allele = 0
                            # templates 1 and 2 cover each of the sites in the 4 SNP long default templates
                            # TODO fix this for other templates (store beginning of missing run rather than the length)
                            if t == 0 or t == 1:
                                missing[haplotype_index] += 1
                        else:
                            # reset the counter of missing sites for this haplotype
                            missing[haplotype_index] = 0

                        last_allele = allele
                        alleles[i] = allele

                        if allele == 0:
                            ppa_0[ppa_0_i] = haplotype_index
                            ppa_0_i += 1
                            div_0[div_0_i] = start_0
                            div_0_i += 1
                            num_matches_end_0 += 1
                            start_0 = 0
                        else:
                            ppa_1[ppa_1_i] = haplotype_index
                            ppa_1_i += 1
                            div_1[div_1_i] = start_1
                            div_1_i += 1
                            num_matches_end_1 += 1
                            start_1 = 0

                    # for matches within the last block of haplotypes
                    if compute_ibd:
                        end = last_site
                        for u in xrange(i - (num_matches_end_0 + num_matches_end_1), i):
                            start = 0
                            for y in xrange(u + 1, i + 1):
                                if div[t][y] > start:
                                    start = div[t][y]
                                match_long_enough = end > start and end - start >= L_m
                                if match_long_enough:
                                    if alleles[u] != alleles[y]:
                                        out_of_sample_match = False
                                        if not in_sample_analysis:
                                            out_of_sample_match = (((ppa[t][u] < num_all_x_haplotypes) and \
                                                                    (ppa[t][y] >= num_all_x_haplotypes)) or \
                                                                   ((ppa[t][y] < num_all_x_haplotypes) and \
                                                                    (ppa[t][u] >= num_all_x_haplotypes)))
                                        if in_sample_analysis or out_of_sample_match:
                                            self.report_templated_segment(ppa[t][u], ppa[t][y], 
                                                    start, end, L_m, L_f, use_phase_correction)
                                else:
                                    break
             
                    # new arrays for next iteration of this template
                    if ppa_0_i > 0:
                        memcpy(&(ppa[t][0]), &ppa_0[0], sizeof(long) * ppa_0_i)
                    if ppa_1_i > 0:
                        memcpy(&(ppa[t][ppa_0_i]), &ppa_1[0], sizeof(long) * ppa_1_i)
                    if div_0_i > 0:
                        memcpy(&(div[t][0]), &div_0[0], sizeof(uint32_t) * div_0_i)
                    if div_1_i > 0:
                        memcpy(&(div[t][div_0_i]), &div_1[0], sizeof(uint32_t) * div_1_i)
            
            if current_site < num_sites_x - 1:
                template_pos += 1
                if template_pos == self.template_length:
                    template_pos = 0

        if compute_ibd: 
            # for trailing matches that extend the length of the alignment
            # iterate over templates
            for t in xrange(self.num_templates):

                # determine which templates to process for this site:
                if self.template_ptr[t][template_pos] == 1:
                    for i in xrange(total_num_active_haplotypes - 1):
                        index_0 = ppa[t][i]
                        start = 0
                        for j in xrange(i + 1, total_num_active_haplotypes):
                            index_1 = ppa[t][j]
                            if div[t][j] > start:
                                start = div[t][j]
                            match_long_enough = current_site > start and current_site - start >= L_m
                            if match_long_enough:
                                out_of_sample_match = False
                                if not in_sample_analysis:
                                    out_of_sample_match = (((index_0 < num_all_x_haplotypes) and \
                                                            (index_1 >= num_all_x_haplotypes)) or \
                                                           ((index_1 < num_all_x_haplotypes) and \
                                                            (index_0 >= num_all_x_haplotypes)))
                                if in_sample_analysis or out_of_sample_match:
                                    self.report_templated_segment(index_0, index_1, start, 
                                            current_site, L_m, L_f, use_phase_correction)
                            else:
                                break
            
            # finally report any last segments
            self.flush_templated_segments(L_f)
       
        if compressed_out_path != NULL:
            fclose(f_out)
        if segments_out_path != NULL:
            fclose(self.segment_output_file)
            self.segment_output_file = NULL

        # free memory
        for i in xrange(self.num_templates):
            free(ppa[i])
            free(div[i])
        free(ppa)
        free(ppa_0)
        free(ppa_1)
        free(div)
        free(div_0)
        free(div_1)
        free(missing)
        free(alleles)
        if vcf_line:
            free(vcf_line)
    
    
    cdef uint8_t report_templated_segment(self, long id1, long id2, uint32_t start, uint32_t end, uint32_t L_m, double L_f, bint use_phase_correction):
        """
        Merges and extends segments found by under each template. When use_phase_correction == True
        a heuristic is used to extend segments across switch errors.
        """
        cdef uint32_t gap_threshold = self.template_length
        gap_threshold += L_m

        cdef long id_temp
        cdef uint8_t out_of_sample = 0
        if id2 < id1:
            id_temp = id1
            id1 = id2
            id2 = id_temp
        if id2 >= self.num_all_x_haplotypes:
            out_of_sample = 1
            id2 = id2 - self.num_all_x_haplotypes

        cdef long id1_com
        cdef long id2_com

        # phase correction heuristic
        if use_phase_correction:


            # get complement haplotypes from the two individuals
            id1_com = self.x_haplotypes.get_haplotype_complement(id1)
            if out_of_sample == 0:
                id2_com = self.x_haplotypes.get_haplotype_complement(id2)
            else:
                id2_com = self.y_haplotypes.get_haplotype_complement(id2)
       
            self_ibd = (id2 == id1_com)
            
            # if either/both individuals are in a switch error then swap haplotypes
            if self.current_switch_error[id1] == 1:
                id_temp = id1
                id1 = id1_com
                id1_com = id_temp
            if (out_of_sample == 0 and self.current_switch_error[id2] == 1) or \
               (out_of_sample == 1 and \
                self.current_switch_error[id2 + self.num_all_x_haplotypes] == 1):
                id_temp = id2
                id2 = id2_com
                id2_com = id_temp
            
            # only check for phase switch if we can't extend a current match
            found_switch_in_id1 = False
            found_switch_in_id2 = False
            if start > self.current_matches_end[id1][id2] + gap_threshold: 
            
                # first check for phase switch by checking if the new segment starts
                # near the end of an old segment on a complement haplotype

                # check if this segment begins at a putative switch site in both individuals:
                if not self_ibd \
                   and start > self.current_matches_end[id1][id2_com] + gap_threshold \
                   and start > self.current_matches_end[id1_com][id2] + gap_threshold \
                   and start <= self.current_matches_end[id1_com][id2_com] + gap_threshold \
                   and start > self.current_matches_end[id1_com][id2_com] - gap_threshold:
                    found_switch_in_id1 = True
                    found_switch_in_id2 = True

                # check if this segment begins at a putative switch site in individual 1:
                elif not self_ibd \
                   and start > self.current_matches_end[id1_com][id2_com] + gap_threshold \
                   and start > self.current_matches_end[id1][id2_com] + gap_threshold \
                   and start <= self.current_matches_end[id1_com][id2] + gap_threshold \
                   and start > self.current_matches_end[id1_com][id2] - gap_threshold:
                    found_switch_in_id1 = True

                # check if this segment begins at a putative switch site in individual 2:
                elif not self_ibd \
                   and start > self.current_matches_end[id1_com][id2_com] + gap_threshold \
                   and start > self.current_matches_end[id1_com][id2] + gap_threshold \
                   and start <= self.current_matches_end[id1][id2_com] + gap_threshold \
                   and start > self.current_matches_end[id1][id2_com] - gap_threshold:
                    found_switch_in_id2 = True
            
                # finally also check for phase switch errors by checking if the new segment is
                # invalidly phased: an invalidly phased segment overlaps with another segment 
                # on the complement of one haplotype but not the other.
                elif start < self.current_matches_end[id1_com][id2]:
                    found_switch_in_id1 = True
                elif start < self.current_matches_end[id1][id2_com]:
                    found_switch_in_id2 = True
            
                if found_switch_in_id1:
                    # swap ids
                    id_temp = id1
                    id1 = id1_com
                    id1_com = id_temp
                    # change the phase switch indicator
                    if self.current_switch_error[id1] == 0:
                        self.current_switch_error[id1] = 1
                        self.current_switch_error[id1_com] = 1
                    else:
                        self.current_switch_error[id1] = 0
                        self.current_switch_error[id1_com] = 0
                if found_switch_in_id2:
                    # swap ids
                    id_temp = id2
                    id2 = id2_com
                    id2_com = id_temp
                    # change the phase switch indicator
                    if out_of_sample == 0:
                        if self.current_switch_error[id2] == 0:
                            self.current_switch_error[id2] = 1
                            self.current_switch_error[id2_com] = 1
                        else:
                            self.current_switch_error[id2] = 0
                            self.current_switch_error[id2_com] = 0
                    else:
                        if self.current_switch_error[id2 + self.num_all_x_haplotypes] == 0:
                            self.current_switch_error[id2 + self.num_all_x_haplotypes] = 1
                            self.current_switch_error[id2_com + self.num_all_x_haplotypes] = 1
                        else:
                            self.current_switch_error[id2 + self.num_all_x_haplotypes] = 0
                            self.current_switch_error[id2_com + self.num_all_x_haplotypes] = 0
            
            # truncate any leftover invalidly phased segments
            if start < self.current_matches_end[id1][id2_com]:
                if start < self.current_matches_start[id1][id2_com]:
                    self.current_matches_start[id1][id2_com] = 0
                    self.current_matches_end[id1][id2_com] = 0
                else:
                    self.current_matches_end[id1][id2_com] = start
            if start < self.current_matches_end[id1_com][id2]:
                if start < self.current_matches_start[id1_com][id2]:
                    self.current_matches_start[id1_com][id2] = 0
                    self.current_matches_end[id1_com][id2] = 0
                else:
                    self.current_matches_end[id1_com][id2] = start
          
        if self.current_matches_end[id1][id2] == 0:
            # the first segment for this pair, so save it
            self.current_matches_start[id1][id2] = start
            self.current_matches_end[id1][id2] = end
        elif start > self.current_matches_end[id1][id2] + gap_threshold:
            # a new segment starts past the gap_threshold, so report the old segment
            self.report_segment(id1, id2, self.current_matches_start[id1][id2], \
                                self.current_matches_end[id1][id2], L_f) 
            self.current_matches_start[id1][id2] = start
            self.current_matches_end[id1][id2] = end
        else:
            # merge the segments
            if start < self.current_matches_start[id1][id2]:
                self.current_matches_start[id1][id2] = start
            if end > self.current_matches_end[id1][id2]:
                self.current_matches_end[id1][id2] = end

        return(0)

    
    cdef uint8_t flush_templated_segments(self, double L_f):
        
        cdef uint32_t id1
        cdef uint32_t id2
        if self.num_all_y_haplotypes == 0:
            for id1 in range(self.num_all_x_haplotypes - 1):
                for id2 in range(id1 + 1, self.num_all_x_haplotypes):
                    if self.current_matches_end[id1][id2] != 0:
                        self.report_segment(id1, id2, self.current_matches_start[id1][id2], \
                                            self.current_matches_end[id1][id2], L_f) 
                        self.current_matches_start[id1][id2] = 0
                        self.current_matches_end[id1][id2] = 0
        else:
            for id1 in range(self.num_all_x_haplotypes):
                for id2 in range(self.num_all_y_haplotypes):
                    if self.current_matches_end[id1][id2] != 0:
                        self.report_segment(id1, id2, self.current_matches_start[id1][id2], \
                                            self.current_matches_end[id1][id2], L_f) 
                        self.current_matches_start[id1][id2] = 0
                        self.current_matches_end[id1][id2] = 0
        return(0)


    cdef uint8_t report_segment(self, long id1, long id2, uint32_t start, uint32_t end, double L_f=3.0):
        """
        Reports segments to output. For in-sample analyses id1 and id2 are sorted
        lexicographically. For out-of-sample analyses id1 is from haplotype alignment X
        and id2 is from haplotype alignment Y.
        """
        cdef double start_cm = self.x_haplotypes.get_genetic_position_from_site(start)
        cdef double end_cm = self.x_haplotypes.get_genetic_position_from_site(end)
        cdef uint32_t start_bp = self.x_haplotypes.get_physical_position_from_site(start)
        cdef uint32_t end_bp = self.x_haplotypes.get_physical_position_from_site(end)
        cdef long id1_name
        cdef long id2_name
        cdef uint8_t id1_hap
        cdef uint8_t id2_hap
        cdef uint8_t out_of_sample = 0

        if end_cm - start_cm >= L_f:

            if self.num_all_y_haplotypes > 0:
                out_of_sample = 1
                id2_name = self.y_haplotypes.get_haplotype_name(id2)
                id2_hap = self.y_haplotypes.get_haplotype_number(id2)
            else:
                id2_name = self.x_haplotypes.get_haplotype_name(id2)
                id2_hap = self.x_haplotypes.get_haplotype_number(id2)
            id1_name = self.x_haplotypes.get_haplotype_name(id1)
            id1_hap = self.x_haplotypes.get_haplotype_number(id1)

            if self.segment_output_file != NULL:
                # write segments to file
                row = ''
                if out_of_sample or id1_name < id2_name:
                    row += str(id1_name) + ','
                    row += str(id2_name) + ','
                    row += str(id1_hap) + ','
                    row += str(id2_hap) + ','
                else:
                    row += str(id2_name) + ','
                    row += str(id1_name) + ','
                    row += str(id2_hap) + ','
                    row += str(id1_hap) + ','
                row += str(start) + ','
                row += str(end) + ','
                #row += str(round(start_cm, 2)) + ','
                #row += str(round(end_cm, 2)) + ','
                row += str(start_cm) + ','
                row += str(end_cm) + ','
                row += str(start_bp) + ','
                row += str(end_bp) + ','
                row += self.chromosome + '\n'
                row_len = len(row)
                temp_row = row.encode()
                fwrite(<char*>temp_row, sizeof(char), row_len, self.segment_output_file)
                self.num_segments += 1
            else:
                # keep segments in memory 
                if out_of_sample or id1_name < id2_name:
                    self.segments_id1[self.num_segments] = id1_name
                    self.segments_id2[self.num_segments] = id2_name
                    self.segments_id1_haplotype[self.num_segments] = id1_hap
                    self.segments_id2_haplotype[self.num_segments] = id2_hap
                else:
                    self.segments_id1[self.num_segments] = id2_name
                    self.segments_id2[self.num_segments] = id1_name
                    self.segments_id1_haplotype[self.num_segments] = id2_hap
                    self.segments_id2_haplotype[self.num_segments] = id1_hap
                self.segments_start[self.num_segments] = start
                self.segments_end[self.num_segments] = end
                self.segments_start_cm[self.num_segments] = start_cm
                self.segments_end_cm[self.num_segments] = end_cm
                self.segments_start_bp[self.num_segments] = start_bp
                self.segments_end_bp[self.num_segments] = end_bp
                self.segments_chromosome.append(self.chromosome)
                self.num_segments += 1
                if self.num_segments == self.segments_capacity:
                    print('Increasing size of segment vectors...')
                    sys.stdout.flush()
                    self.segments_capacity *= 2
                    self.segments_id1 = <long *>realloc(self.segments_id1, self.segments_capacity * sizeof(long))
                    self.segments_id2 = <long *>realloc(self.segments_id2, self.segments_capacity * sizeof(long))
                    self.segments_id1_haplotype = <uint8_t *>realloc(self.segments_id1_haplotype, self.segments_capacity * sizeof(uint8_t))
                    self.segments_id2_haplotype = <uint8_t *>realloc(self.segments_id2_haplotype, self.segments_capacity * sizeof(uint8_t))
                    self.segments_start = <uint32_t *>realloc(self.segments_start, self.segments_capacity * sizeof(uint32_t))
                    self.segments_end = <uint32_t *>realloc(self.segments_end, self.segments_capacity * sizeof(uint32_t))
                    self.segments_start_cm = <double *>realloc(self.segments_start_cm, self.segments_capacity * sizeof(double))
                    self.segments_end_cm = <double *>realloc(self.segments_end_cm, self.segments_capacity * sizeof(double))
                    self.segments_start_bp = <uint32_t *>realloc(self.segments_start_bp, self.segments_capacity * sizeof(uint32_t))
                    self.segments_end_bp = <uint32_t *>realloc(self.segments_end_bp, self.segments_capacity * sizeof(uint32_t))
        return(0)


    cdef get_segments(self):
        cdef uint64_t i = 0
        id1 = []
        id2 = []
        id1_haplotype = []
        id2_haplotype = []
        start = []
        end = []
        start_cm = []
        end_cm = []
        start_bp = []
        end_bp = []
        for i in xrange(self.num_segments):
            id1.append(self.segments_id1[i])
            id2.append(self.segments_id2[i])
            id1_haplotype.append(self.segments_id1_haplotype[i])
            id2_haplotype.append(self.segments_id2_haplotype[i])
            start.append(self.segments_start[i])
            end.append(self.segments_end[i])
            start_cm.append(self.segments_start_cm[i])
            end_cm.append(self.segments_end_cm[i])
            start_bp.append(self.segments_start_bp[i])
            end_bp.append(self.segments_end_bp[i])
        chromosome = np.array(self.segments_chromosome, dtype=object)
        id1 = np.array(id1, dtype=np.int64)
        id2 = np.array(id2, dtype=np.int64)
        id1_haplotype = np.array(id1_haplotype, dtype=np.int64)
        id2_haplotype = np.array(id2_haplotype, dtype=np.int64)
        start = np.array(start, dtype=np.int64)
        end = np.array(end, dtype=np.int64)
        start_cm = np.array(start_cm, dtype=np.float64)
        end_cm = np.array(end_cm, dtype=np.float64)
        start_bp = np.array(start_bp, dtype=np.int64)
        end_bp = np.array(end_bp, dtype=np.int64)
        self.ibd_results = pd.DataFrame({'chromosome': chromosome,
                                         'id1': id1,
                                         'id2': id2,
                                         'id1_haplotype': id1_haplotype,
                                         'id2_haplotype': id2_haplotype,
                                         'start': start,
                                         'end': end,
                                         'start_cm': start_cm,
                                         'end_cm': end_cm,
                                         'start_bp': start_bp,
                                         'end_bp': end_bp})
        return self.ibd_results
    
    
    cdef clear_segment_lists(self):
        # make all matching segment lists empty
        self.num_segments = 0
        self.segments_capacity = 10000000
        free(self.segments_id1)
        free(self.segments_id2)
        free(self.segments_id1_haplotype)
        free(self.segments_id2_haplotype)
        free(self.segments_start)
        free(self.segments_end)
        free(self.segments_start_cm)
        free(self.segments_end_cm)
        free(self.segments_start_bp)
        free(self.segments_end_bp)
        self.segments_id1 = <long *>malloc(self.segments_capacity * sizeof(long))
        self.segments_id2 = <long *>malloc(self.segments_capacity * sizeof(long))
        self.segments_id1_haplotype = <uint8_t *>malloc(self.segments_capacity * sizeof(uint8_t))
        self.segments_id2_haplotype = <uint8_t *>malloc(self.segments_capacity * sizeof(uint8_t))
        self.segments_start = <uint32_t *>malloc(self.segments_capacity * sizeof(uint32_t))
        self.segments_end = <uint32_t *>malloc(self.segments_capacity * sizeof(uint32_t))
        self.segments_start_cm = <double *>malloc(self.segments_capacity * sizeof(double))
        self.segments_end_cm = <double *>malloc(self.segments_capacity * sizeof(double))
        self.segments_start_bp = <uint32_t *>malloc(self.segments_capacity * sizeof(uint32_t))
        self.segments_end_bp = <uint32_t *>malloc(self.segments_capacity * sizeof(uint32_t))
        self.segments_chromosome = []


    
   

cdef class TPBWTAnalysis:
    """
    Class that manages TPBWT analyses over multiple chromosomes.
    Setting template to [[1]] collapses the TPBWT down to the PBWT.
    """

    def __init__(self, template=None):

        self.tpbwt = TPBWT(template) 
        self.x_haplotypes = None
        self.y_haplotypes = None


    cpdef compress_alignment(self, compressed_out_path, HaplotypeAlignment x_haplotypes, HaplotypeAlignment y_haplotypes=None, bint verbose=True):
        """
        Method that TPBWT-compresses a haplotype alignment without simultaneously computing IBD
        so it uses less memory. If two alignments are input then it compresses them into a single
        combined haplotype alignment.
        """
        start_time = timeit.default_timer()
        self.x_haplotypes = x_haplotypes
        self.y_haplotypes = y_haplotypes
        if not os.path.isdir(compressed_out_path):
            os.mkdir(compressed_out_path)
        cdef uint8_t c
        chromosomes = self.x_haplotypes.get_chromosomes()
        for c in xrange(len(chromosomes)):
            chromo = chromosomes[c]
            if verbose:
                print('\nTPBWT-compressing chromosome ' + chromo + '...')
            sys.stdout.flush()
            self.x_haplotypes.set_chromosome(chromo)
            if self.y_haplotypes is not None:
                self.y_haplotypes.set_chromosome(chromo)
            chromo_out_path = os.path.join(compressed_out_path, chromo + '.tpbwt').encode()
            self.tpbwt.transform(x_haplotypes, y_haplotypes, chromo_out_path, 300, 3.0, 10, False, NULL, verbose, False)
            if verbose:
                print('Writing haplotype metadata to binary file...')
            sys.stdout.flush()
            out_haplo_path = os.path.join(compressed_out_path, chromo + '_md.tpbwt').encode()
            self.tpbwt.write_haplotype_data(out_haplo_path)
        end_time = timeit.default_timer()
        if verbose:
            print('Done. Time elapsed: ' + str(round(end_time - start_time, 2)))
            sys.stdout.flush()
        return(0)


    cpdef compute_ibd(self, HaplotypeAlignment x_haplotypes, HaplotypeAlignment y_haplotypes=None, compressed_out_path=None, uint32_t L_m=300, double L_f=3.0, uint32_t missing_site_threshold=10, bint use_phase_correction=True, chromosome='all', segments_out_path=None, bint verbose=True):
        """
        Method that configures the output format of the IBD compute and then iterates over each 
        chromosome. For each chromosome it sets the current chromosome in the haplotype 
        alignments and runs TPBWT. Optionally TPBWT-compresses while computing IBD.
        """
        start_time = timeit.default_timer()
        self.x_haplotypes = x_haplotypes
        self.y_haplotypes = y_haplotypes
        if compressed_out_path:
            if not os.path.isdir(compressed_out_path):
                os.mkdir(compressed_out_path)
        if segments_out_path:
            if not os.path.isdir(segments_out_path):
                os.mkdir(segments_out_path)
        cdef uint8_t c
        if chromosome == 'all':
            chromosomes = self.x_haplotypes.get_chromosomes()
        else:
            chromosomes = [chromosome]
        for c in xrange(len(chromosomes)):
            chromo = chromosomes[c]
            if verbose:
                print('\nIBD compute for chromosome ' + chromo + '...')
            sys.stdout.flush()
            self.x_haplotypes.set_chromosome(chromo) 
            if self.y_haplotypes is not None:
                self.y_haplotypes.set_chromosome(chromo)
            if not compressed_out_path and not segments_out_path:
                self.tpbwt.transform(x_haplotypes, y_haplotypes, NULL, L_m, L_f, missing_site_threshold, use_phase_correction, NULL, verbose, True)
            elif compressed_out_path and not segments_out_path:
                chromo_out_path = os.path.join(compressed_out_path, chromo + '.tpbwt').encode()
                self.tpbwt.transform(x_haplotypes, y_haplotypes, chromo_out_path, L_m, L_f, missing_site_threshold, use_phase_correction, NULL, verbose, True)
            elif not compressed_out_path and segments_out_path:
                chromo_out_path = os.path.join(segments_out_path, chromo + '.csv').encode()
                self.tpbwt.transform(x_haplotypes, y_haplotypes, NULL, L_m, L_f, missing_site_threshold, use_phase_correction, chromo_out_path, verbose, True)
            elif compressed_out_path and segments_out_path:
                chromo_out_path1 = os.path.join(compressed_out_path, chromo + '.tpbwt').encode()
                chromo_out_path2 = os.path.join(segments_out_path, chromo + '.csv').encode()
                self.tpbwt.transform(x_haplotypes, y_haplotypes, chromo_out_path1, L_m, L_f, missing_site_threshold, use_phase_correction, chromo_out_path2, verbose, True)
            if compressed_out_path:
                if verbose:
                    print('Writing haplotype metadata to binary file...')
                sys.stdout.flush()
                out_haplo_path = os.path.join(compressed_out_path, chromo + '_md.tpbwt').encode()
                self.tpbwt.write_haplotype_data(out_haplo_path)
        if verbose:
            print('\nDone computing IBD segments. Finishing up...')
            print('Number of IBD segments found = ' + str(self.tpbwt.num_segments))
            sys.stdout.flush()
        if not segments_out_path:
            if verbose:
                print('Building final dataframe...')
                sys.stdout.flush()
            results = self.tpbwt.get_segments()
        else:
            results = 0
        self.tpbwt.clear_segment_lists()
        end_time = timeit.default_timer()
        if verbose:
            print('Done. Time elapsed: ' + str(round(end_time - start_time, 2)))
            sys.stdout.flush()
        return(results)
