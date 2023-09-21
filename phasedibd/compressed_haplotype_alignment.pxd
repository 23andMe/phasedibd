from libc.stdint cimport uint8_t, uint16_t, uint32_t
from libc.stdio cimport *
from .haplotype_alignment cimport HaplotypeAlignment


cdef class CompressedHaplotypeAlignment(HaplotypeAlignment):

    cdef double [:] chromo_genetic_position
    cdef uint32_t [:] chromo_physical_position
    cdef path
    cdef file_extension
    cdef long [:] haplotype_names_view
    cdef uint8_t [:] haplotype_numbers_view
    cdef FILE* compressed_tpbwt_file
    cdef uint8_t current_allele
    cdef uint16_t current_allele_total
    cdef uint16_t current_allele_count
    cdef uint32_t haplotype_count
      
    cdef double get_genetic_position_from_site(self, long site)
    cdef uint32_t get_physical_position_from_site(self, long site)
    cdef get_chromosomes(self)
    cdef get_haplotype_sex(self, long haplotype_index)
    cdef uint32_t get_num_sites(self, chromosome=?)
    cdef long get_haplotype_index(self, long haplotype_name)
    cdef long get_haplotype_name(self, long haplotype_index)
    cdef uint8_t get_allele(self, uint32_t haplotype_index, uint32_t alignment_index)
    cdef uint8_t get_haplotype_number(self, long haplotype_index)
    cdef uint8_t read_haplotype_data(self, path)
    cdef set_chromosome(self, chromosome)
    cdef uint8_t** get_alignment(self)
    cdef FILE* get_tpbwt_file(self)
