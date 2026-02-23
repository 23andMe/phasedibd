
from libc.stdio cimport *
from phasedibd.haplotype_alignment cimport HaplotypeAlignment
from phasedibd.compressed_haplotype_alignment cimport CompressedHaplotypeAlignment
from phasedibd.vcf_haplotype_alignment cimport VcfHaplotypeAlignment
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t


cdef class TPBWT:

    cdef HaplotypeAlignment x_haplotypes
    cdef HaplotypeAlignment y_haplotypes
    cdef template
    cdef uint8_t** template_ptr
    cdef uint8_t num_templates
    cdef uint8_t template_length
    cdef public chromosome
    cdef uint32_t num_all_x_haplotypes
    cdef uint32_t num_all_y_haplotypes
    cdef long* segments_id1
    cdef long* segments_id2
    cdef uint8_t* segments_id1_haplotype
    cdef uint8_t* segments_id2_haplotype
    cdef uint32_t* segments_start
    cdef uint32_t* segments_end
    cdef double* segments_start_cm
    cdef double* segments_end_cm
    cdef uint32_t* segments_start_bp
    cdef uint32_t* segments_end_bp
    cdef uint64_t num_segments
    cdef uint64_t segments_capacity
    cdef public segments_chromosome
    cdef public ibd_results
    cdef uint32_t** current_matches_start
    cdef uint32_t** current_matches_end
    cdef uint8_t* current_switch_error 
    cdef FILE* segment_output_file
    
    cdef uint8_t allocate_current_match_arrays(self)
    cdef transform(self, HaplotypeAlignment x_haplotypes, HaplotypeAlignment y_haplotypes=?, char* compressed_out_path=?, uint32_t L_m=?, double L_f=?, uint32_t missing_site_threshold=?, bint use_phase_correction=?, char* segments_out_path=?, bint verbose=?, bint compute_ibd=?)
    cdef uint8_t flush_templated_segments(self, double L_f)
    cdef get_segments(self)
    cdef clear_segment_lists(self)
    cdef uint8_t report_segment(self, long id1, long id2, uint32_t start, uint32_t end, double L_f=?)
    cdef uint8_t report_templated_segment(self, long id1, long id2, uint32_t start, uint32_t end, uint32_t L_m, double L_f, bint use_phase_correction)
    cdef uint8_t write_haplotype_data(self, path)


cdef class TPBWTAnalysis:

    cdef HaplotypeAlignment x_haplotypes
    cdef HaplotypeAlignment y_haplotypes
    cdef TPBWT tpbwt

    cpdef compress_alignment(self, compressed_out_path, HaplotypeAlignment x_haplotypes, HaplotypeAlignment y_haplotypes=?, bint verbose=?)
    cpdef compute_ibd(self, HaplotypeAlignment x_haplotypes, HaplotypeAlignment y_haplotypes=?, compressed_out_path=?, uint32_t L_m=?, double L_f=?, uint32_t missing_site_threshold=?, bint use_phase_correction=?, chromosome=?, segments_out_path=?, bint verbose=?)

