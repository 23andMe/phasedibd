# cython: profile=True

from libc.stdlib cimport malloc, free
from haplotype_alignment cimport HaplotypeAlignment
import numpy as np
cimport numpy as np
import pandas as pd




cdef class VcfHaplotypeAlignment(HaplotypeAlignment):


    def __init__(self, vcf_file, map_file=None):

        self.current_chromosome = ''
        self.vcf_file = vcf_file
        self.num_all_haplotypes = 0
        self.num_active_haplotypes = 0
        self.num_sites = 0
        self.header_length = 0
        self.map_file = map_file
        self.physical_positions = []
        self.load_data()    
        self.current_site = -1
  

    cdef void append_physical_positions(self, int pos):
        self.physical_positions.append(pos)
        return


    cdef size_t get_header_length(self):
        return(self.header_length)


    cdef size_t get_sample_start(self):
        return(self.sample_start)


    cdef char* get_vcf_file(self):
        return(self.vcf_file.encode())


    cdef uint8_t load_data(self):

        sample_names = []
        done_reading_header = False
        setting_sample_names = False

        # scan through vcf to get sample names and num sites
        f = fopen(self.vcf_file.encode(), 'r')
        cdef char *line = NULL
        cdef char *word
        cdef size_t n = 0
        cdef ssize_t line_len = getline(&line, &n, f)
        cdef size_t word_num = 0
        cdef size_t i = 0
        self.sample_start = 0
        vcf_delimiter_py = '\t'.encode()
        cdef char *vcf_delimiter = vcf_delimiter_py
        while line_len != -1: 
            if not done_reading_header:
                # parse line
                self.header_length += 1
                word = strtok(line, vcf_delimiter)
                word_num = 0
                while word != NULL:
                    if word.decode() == 'FORMAT':
                        setting_sample_names = True
                        self.sample_start = word_num + 1
                        current_sample = 0
                    elif setting_sample_names:
                        # get sample names
                        #sample_names.append(word.decode().rstrip())
                        sample_names.append(current_sample)
                        current_sample += 1
                    word = strtok(NULL, vcf_delimiter)
                    word_num += 1
            else:
                self.num_sites += 1
                # the first word is the chromosome
                word = strtok(line, vcf_delimiter)
                if self.num_sites == 1:
                    self.current_chromosome = word.decode().rstrip()
                # second is the physical position
                word = strtok(NULL, vcf_delimiter)
                self.physical_positions.append(int(word.decode().rstrip()))

            if setting_sample_names:
                # allocate haplotype metadata arrays; currently assumes diploid
                self.num_all_haplotypes = 2 * len(sample_names)
                self.num_active_haplotypes = self.num_all_haplotypes
                self.haplotype_names_view = np.empty(shape=(self.num_all_haplotypes), dtype=np.int64)
                self.haplotype_numbers_view = np.empty(shape=(self.num_all_haplotypes), dtype=np.uint8)
                self.haplotypes_pointer = <uint8_t **>malloc(self.num_all_haplotypes * sizeof(uint8_t*))
                self.active_haplotype_indices = <long *>malloc(self.num_active_haplotypes * sizeof(long))

                for i in xrange(len(sample_names)):
                    name = int(sample_names[i])
                    self.haplotype_names_view[i * 2] = name
                    self.haplotype_names_view[(i * 2) + 1] = name
                    self.haplotype_numbers_view[i * 2] = 0
                    self.haplotype_numbers_view[(i * 2) + 1] = 1
                    self.active_haplotype_indices[i * 2] = i * 2
                    self.active_haplotype_indices[(i * 2) + 1] = (i * 2) + 1

                    # allocate and fill vectors to hold haplotype data for a single site
                    self.haplotypes_pointer[i * 2] = <uint8_t *>malloc(1 * sizeof(uint8_t))
                    self.haplotypes_pointer[(i * 2) + 1] = <uint8_t *>malloc(1 * sizeof(uint8_t))
                    self.haplotypes_pointer[i * 2][0] = 0
                    self.haplotypes_pointer[(i * 2) + 1][0] = 0
                done_reading_header = True
                setting_sample_names = False 

            # read next line
            line_len = getline(&line, &n, f)
      
        fclose(f)
        if line:
            free(line)
        
        self.chromo_physical_position = np.asarray(self.physical_positions, dtype=np.uint32)

        # set genetic map
        if self.map_file is None:
            # set genetic positions 
            self.chromo_genetic_position = np.asarray(self.physical_positions)/1000000.0
        else:
            # read PLINK format genetic map file setting physical and genetic positions
            m = pd.read_csv(self.map_file, header=None, names=['chromo','id','cm','bp'], sep=r'\s+')
            if np.array(m['bp']).shape != np.array(self.physical_positions).shape or \
                    (np.array(m['bp']) != np.array(self.physical_positions)).any():
                raise Exception('The sites in the VCF do not match the sites in the genetic map.')
            self.chromo_genetic_position = np.array(m['cm'])


    cdef uint8_t** get_alignment(self):
        return(self.haplotypes_pointer)
   
    
    cdef uint8_t get_allele(self, uint32_t haplotype_index, uint32_t alignment_index):
        return(self.haplotypes_pointer[haplotype_index][alignment_index])

    
    cdef get_chromosomes(self):
        return([self.current_chromosome])


    cdef set_chromosome(self, chromosome):
        return()


    cdef long* get_active_haplotype_indices(self):
        return(self.active_haplotype_indices)


    cdef double get_genetic_position_from_site(self, long site):
        return(self.chromo_genetic_position[site])
    
    
    cdef uint32_t get_physical_position_from_site(self, long site):
        return(self.chromo_physical_position[site])
    
    
    cdef long get_haplotype_index(self, long haplotype_name):
        cdef long i
        for i in xrange(len(self.haplotype_names_view)):
            if haplotype_name == self.haplotype_names_view[i]:
                return(i)


    cdef long get_haplotype_name(self, long haplotype_index):
        """
        Return the id of the individual who has this haplotype.
        """
        return(self.haplotype_names_view[haplotype_index])
    
   
    cdef uint8_t get_haplotype_number(self, long haplotype_index):
        """
        Is this haplotype 0 or 1 for the individual?
        """
        return(self.haplotype_numbers_view[haplotype_index])
    
    
    cdef uint32_t get_num_all_haplotypes(self):
        return(self.num_all_haplotypes)


    cdef void update_position(self, long site):

        cdef size_t num_read_lines = 0
        cdef char *line = NULL
        cdef char *word
        cdef size_t n = 0
        cdef ssize_t line_len = 0
        cdef size_t word_num = 0
        cdef size_t current_sample = 0
        vcf_delimiter_py = '\t'.encode()
        cdef char *vcf_delimiter = vcf_delimiter_py

        if site != self.current_site:
            # read next line of vcf
            if self.current_site == -1:
                self.file = fopen(self.vcf_file.encode(), 'r')
                while num_read_lines <= self.header_length:
                    line_len = getline(&line, &n, self.file)
                    num_read_lines += 1
            else:
                line_len = getline(&line, &n, self.file)
            self.current_site = site

            if line_len != -1: 
                # parse line
                word = strtok(line, vcf_delimiter)
                while word != NULL:
                    if word_num >= self.sample_start:
                        # get haplotype alleles for this site
                        # 48 is ascii for '0'
                        # 49 is ascii for '1'
                        if word[0] == 48:
                            self.haplotypes_pointer[current_sample * 2][0] = 0
                        elif word[0] == 49:
                            self.haplotypes_pointer[current_sample * 2][0] = 1
                        else:
                            self.haplotypes_pointer[current_sample * 2][0] = 2
                        if word[2] == 48:
                            self.haplotypes_pointer[(current_sample * 2) + 1][0] = 0
                        elif word[2] == 49:
                            self.haplotypes_pointer[(current_sample * 2) + 1][0] = 1
                        else:
                            self.haplotypes_pointer[(current_sample * 2) + 1][0] = 2
                        current_sample += 1
                    word = strtok(NULL, vcf_delimiter)
                    word_num += 1
            else:
                fclose(self.file)
            if line:
                free(line)
            
