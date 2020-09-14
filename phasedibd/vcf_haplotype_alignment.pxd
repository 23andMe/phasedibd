from libc.stdint cimport uint8_t, uint32_t
from libc.stdio cimport *
from libc.stdlib cimport *
from libc.string cimport *
from haplotype_alignment cimport HaplotypeAlignment


cdef class VcfHaplotypeAlignment(HaplotypeAlignment):

    cdef vcf_file
    cdef map_file
    cdef long [:] haplotype_names_view
    cdef uint8_t [:] haplotype_numbers_view
    cdef double [:] chromo_genetic_position
    cdef uint32_t [:] chromo_physical_position
    cdef physical_positions
    cdef size_t sample_start 
    cdef size_t header_length
    cdef FILE* file 
    
    cdef void append_physical_positions(self, int pos)
    cdef size_t get_header_length(self)
    cdef size_t get_sample_start(self)
    cdef char* get_vcf_file(self)
    cdef uint8_t load_data(self)
      
    cdef uint8_t get_allele(self, uint32_t haplotype_index, uint32_t alignment_index)
    cdef uint8_t** get_alignment(self)
    cdef get_chromosomes(self)
    
    cdef set_chromosome(self, chromosome)
    cdef long* get_active_haplotype_indices(self)
    cdef double get_genetic_position_from_site(self, long site)
    cdef uint32_t get_physical_position_from_site(self, long site)
    cdef long get_haplotype_index(self, long haplotype_name)
    cdef long get_haplotype_name(self, long haplotype_index)
    cdef uint8_t get_haplotype_number(self, long haplotype_index)
    cdef uint32_t get_num_all_haplotypes(self)
    cdef void update_position(self, long site)
