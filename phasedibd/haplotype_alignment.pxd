
from libc.stdint cimport uint8_t, uint32_t

cdef class HaplotypeAlignment:

    cdef chromosomes
    cdef current_chromosome
    cdef current_site
    cdef long* active_haplotype_indices
    cdef haplotype_names
    cdef haplotype_numbers
    cdef haplotype_sexes
    cdef haplotypes

    cdef void append_physical_positions(self, int pos)
    cdef size_t get_header_length(self)
    cdef size_t get_sample_start(self)
    cdef char* get_vcf_file(self)
    cdef uint32_t num_all_haplotypes
    cdef uint32_t num_active_haplotypes
    cdef uint32_t num_sites
    cdef uint8_t [:,:] haplotypes_view
    cdef uint8_t** haplotypes_pointer

    cdef double get_genetic_position_from_site(self, long site)
    cdef uint32_t get_physical_position_from_site(self, long site)
    cdef get_chromosomes(self)
    cdef get_current_chromosome(self)
    cdef get_haplotype_sex(self, long haplotype_index)
    cdef uint32_t get_num_all_haplotypes(self)
    cdef uint32_t get_num_active_haplotypes(self)
    cdef uint32_t get_num_sites(self, chromosome=?)
    cdef long* get_active_haplotype_indices(self)
    cdef long get_haplotype_index(self, long haplotype_name)
    cdef long get_haplotype_name(self, long haplotype_index)
    cdef long get_haplotype_complement(self, long haplotype_index)
    cdef uint8_t get_allele(self, uint32_t haplotype_index, uint32_t alignment_index)
    cdef uint8_t get_haplotype_number(self, long haplotype_index)
    cdef set_chromosome(self, chromosome)
    cdef uint8_t** get_alignment(self)
    cdef set_active_haplotype_indices(self)
    cdef void update_position(self, long site)

