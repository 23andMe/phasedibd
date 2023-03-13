# cython: profile=True

import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc, free, realloc

cdef class HaplotypeAlignment:
    """
    Base class for VCFHaplotypeAlignment and FblockHaplotypeAlignment.
    Contains test data.
    """

    def __init__(self, haplotype_array=None, chromosomes=None):
         
        if haplotype_array is None:
            # use test data 
            self.haplotypes = [[0, 1, 0, 1, 0, 1],
                            [1, 1, 2, 0, 0, 1],  # 2 represents missing/invalid allele
                            [1, 1, 1, 1, 1, 1],
                            [0, 1, 1, 1, 1, 0],
                            [0, 0, 0, 0, 0, 0],
                            [1, 0, 0, 0, 1, 0],
                            [1, 1, 0, 0, 0, 1],
                            [0, 1, 0, 1, 1, 0]]
        else:
            self.haplotypes = haplotype_array
        
        if chromosomes is None:
            self.chromosomes = ['1']
        else:
            self.chromosomes = chromosomes

        self.current_chromosome = ''
        
        self.num_all_haplotypes = len(self.haplotypes)
        self.num_active_haplotypes = len(self.haplotypes)
        self.num_sites = len(self.haplotypes[0])
        self.haplotype_names = range(self.num_all_haplotypes)
        self.haplotype_sexes = ['F'] * self.num_all_haplotypes
        self.haplotypes_view = np.array(self.haplotypes, dtype=np.uint8)
        self.haplotypes_pointer = <uint8_t **>malloc(self.num_all_haplotypes * sizeof(uint8_t*))
        for i in xrange(self.num_all_haplotypes):
            self.haplotypes_pointer[i] = &self.haplotypes_view[i][0]


    def __dealloc__(self):
        free(self.haplotypes_pointer)
        free(self.active_haplotype_indices)


    cdef uint8_t** get_alignment(self):
        return(&self.haplotypes_pointer[0])


    cdef uint8_t get_allele(self, uint32_t haplotype_index, uint32_t alignment_index):
        return(self.haplotypes_pointer[haplotype_index][alignment_index])


    cdef get_current_chromosome(self):
        return(self.current_chromosome)


    cdef get_chromosomes(self):
        return(self.chromosomes)
    
    
    cdef double get_genetic_position_from_site(self, long site):
        """
        This function assume the current chromosome was already set by set_chromosome().
        """
        return(site)
    
    
    cdef uint32_t get_physical_position_from_site(self, long site):
        """
        This function assume the current chromosome was already set by set_chromosome().
        """
        return(site)

    
    cdef long* get_active_haplotype_indices(self):
        return(self.active_haplotype_indices)


    cdef long get_haplotype_index(self, long haplotype_name):
        for index, name in self.haplotype_names.iteritems():
            if name == haplotype_name:
                return(index)


    cdef long get_haplotype_name(self, long haplotype_index):
        return(self.haplotype_names[haplotype_index])
    
    
    cdef long get_haplotype_complement(self, long haplotype_index):
        """
        Return the index of the complement haplotype.
        """
        if haplotype_index % 2 == 0:
            return(haplotype_index + 1)
        else:
            return(haplotype_index - 1)
   

    cdef uint8_t get_haplotype_number(self, long haplotype_index):
        """
        Is this haplotype 0 or 1 for the individual?
        """
        return(0)
    
    
    cdef get_haplotype_sex(self, long haplotype_index):
        return(self.haplotype_sexes[haplotype_index])
    

    cdef uint32_t get_num_all_haplotypes(self):
        return(len(self.haplotypes))
    
    
    cdef uint32_t get_num_active_haplotypes(self):
        return(self.num_active_haplotypes)
    
    
    cdef uint32_t get_num_sites(self, chromosome=None):
        return(self.num_sites)

    
    cdef set_chromosome(self, chromosome):
        """
        Set the current chromosome to query.
        """
        if chromosome != self.current_chromosome:
            self.current_chromosome = chromosome
            self.set_active_haplotype_indices()


    cdef set_active_haplotype_indices(self):
        """
        This function determines which haplotypes in the haplotype
        alignment the PBWT should run over. All haplotypes are used
        except in the sex chromosomes, where some haplotypes are not
        included depending on the sex of the individual.
        """
        free(self.active_haplotype_indices)
        cdef long [:] haplos
        cdef int i
        if self.current_chromosome != 'X' and self.current_chromosome != 'Y': 
            self.num_active_haplotypes = self.get_num_all_haplotypes()
            haplos = np.arange(self.num_active_haplotypes, dtype=np.int64)
        else:
            haplotype_indices = []
            if self.current_chromosome == 'X':
                for i in xrange(self.get_num_all_haplotypes()):
                    sex = self.get_haplotype_sex(i)
                    if sex == 'F':
                        haplotype_indices.append(i)
                    elif self.get_haplotype_number(i) == 0:
                        # only add haplotype 0 for males
                        haplotype_indices.append(i)
            elif self.current_chromosome == 'Y':
                for i in xrange(self.get_num_all_haplotypes()):
                    sex = self.get_haplotype_sex(i)
                    if sex == 'M' and self.get_haplotype_number(i) == 0:
                        # only add haplotype 0 for males
                        haplotype_indices.append(i)
            self.num_active_haplotypes = len(haplotype_indices)
            haplos = np.asarray(haplotype_indices, dtype=np.int64)
        self.active_haplotype_indices = <long *>malloc(self.num_active_haplotypes * sizeof(long))
        for i in xrange(self.num_active_haplotypes):
            self.active_haplotype_indices[i] = haplos[i]


    cdef void update_position(self, long site):
        # not implemented; must be cdef-declared for derived classes
        return
    
    cdef void append_physical_positions(self, int pos):
        # not implemented; must be cdef-declared for derived classes
        return
    
    cdef size_t get_header_length(self):
        # not implemented; must be cdef-declared for derived classes
        return(0)
    
    cdef size_t get_sample_start(self):
        # not implemented; must be cdef-declared for derived classes
        return(0)

    cdef char* get_vcf_file(self):
        # not implemented; must be cdef-declared for derived classes
        return('')


