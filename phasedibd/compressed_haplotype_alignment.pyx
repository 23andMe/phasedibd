# cython: profile=True
# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION

from libc.stdlib cimport malloc, free
import numpy as np
cimport numpy as np
from os import listdir
import pandas as pd
from phasedibd.haplotype_alignment cimport HaplotypeAlignment



cdef class CompressedHaplotypeAlignment(HaplotypeAlignment):


    def __init__(self, path):
        """
        path should be a directory containing *_md.tpbwt and *.tpbwt files.
        There should be a _md.tpbwt and .tpbwt file for each chromosome.
        Call set_chromosome() to start loading data.
        """
        self.current_chromosome = ''
        self.path = path
        self.file_extension = '.tpbwt'
   

    cdef uint8_t read_haplotype_data(self, path):
        """
        Reads haplotype names and numbers from binary file.
        """
        # total num haplotypes 
        f = fopen(path.encode(), "rb")
        cdef uint32_t num_all_haplotypes
        fread(&num_all_haplotypes, sizeof(uint32_t), 1, f)
        self.num_all_haplotypes = num_all_haplotypes
         
        # read in names and haplotype #s for all M haplotypes
        self.haplotype_names_view = np.empty(shape=(self.num_all_haplotypes), dtype=np.int64)
        self.haplotype_numbers_view = np.empty(shape=(self.num_all_haplotypes), dtype=np.uint8)
        cdef long temp_name
        cdef uint8_t temp_hap
        for i in xrange(num_all_haplotypes):
            fread(&temp_name, sizeof(long), 1, f)
            fread(&temp_hap, sizeof(uint8_t), 1, f)
            self.haplotype_names_view[i] = temp_name
            self.haplotype_numbers_view[i] = temp_hap
        
        # read num sites, chromosome, genetic/physical positions of sites
        cdef uint32_t num_sites
        fread(&num_sites, sizeof(uint32_t), 1, f)
        self.num_sites = num_sites
        cdef uint8_t chromo_length
        fread(&chromo_length, sizeof(uint8_t), 1, f)
        cdef char* chromo = <char *>malloc(chromo_length * sizeof(char))
        fread(&chromo[0], sizeof(char), chromo_length, f)
        self.current_chromosome = ''
        for i in range(chromo_length):
            self.current_chromosome += chr(chromo[i])
        cdef double gen_pos
        cdef uint32_t phy_pos
        self.chromo_genetic_position = np.empty(shape=(self.num_sites), dtype=np.float64)
        self.chromo_physical_position = np.empty(shape=(self.num_sites), dtype=np.uint32)
        for i in xrange(num_sites):
            fread(&gen_pos, sizeof(double), 1, f)
            self.chromo_genetic_position[i] = gen_pos
            fread(&phy_pos, sizeof(uint32_t), 1, f)
            self.chromo_physical_position[i] = phy_pos
        fclose(f)
            

    cdef uint8_t** get_alignment(self):
        raise NotImplementedError('get_alignment is not implemented for compressed haplotypes.')
    
    
    cdef uint8_t get_allele(self, uint32_t haplotype_index, uint32_t alignment_index):
        """
        This function does not use haplotype_index or alignment_index to look up the
        allele. Instead it assumes the TPBWT-compressed haplotypes were written in the
        correct order to the binary file.
        """
        if self.haplotype_count == 0:
            # start reading in haplotypes for a new column of the alignment
            fread(&self.current_allele, sizeof(uint8_t), 1, self.compressed_tpbwt_file)
            fread(&self.current_allele_total, sizeof(uint16_t), 1, self.compressed_tpbwt_file)
            # keep track of the number of haplotypes in this allele run
            self.current_allele_count = 1
            # keep track of the number of haplotypes seen in this column
            self.haplotype_count = 1
        elif self.current_allele_count == self.current_allele_total:
            # read in the next allele run for this column
            fread(&self.current_allele, sizeof(uint8_t), 1, self.compressed_tpbwt_file)
            fread(&self.current_allele_total, sizeof(uint16_t), 1, self.compressed_tpbwt_file)
            self.current_allele_count = 1
            self.haplotype_count += 1
        else:
            # continue in the current allele run
            self.current_allele_count += 1
            self.haplotype_count += 1
        if self.haplotype_count == self.num_active_haplotypes:
            # we've reached the end of this column of the alignment
            self.haplotype_count = 0

        return(self.current_allele)


    cdef get_chromosomes(self):
        """
        The chromosomes to be analysed are determined by the metadata in the directory.
        """
        metadata_files = [f for f in listdir(self.path) if '_md.' in f]
        chromos = []
        for f in metadata_files:
            pos = f.find('_md')
            chromos.append(f[:pos])
        return(chromos)
    
    
    cdef double get_genetic_position_from_site(self, long site):
        """
        This function assume the current chromosome was already set by set_chromosome().
        Returns the genetic position in cM of a site.
        """
        return(self.chromo_genetic_position[site])
    
    
    cdef uint32_t get_physical_position_from_site(self, long site):
        """
        This function assume the current chromosome was already set by set_chromosome().
        Returns the physical position in bp of a site.
        """
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

    
    cdef get_haplotype_sex(self, long haplotype_index):
        raise NotImplementedError('get_haplotype_sex is not implmented for compressed haplotypes.')

    
    cdef uint32_t get_num_all_haplotypes(self):
        return(self.num_all_haplotypes)
    
    
    cdef uint32_t get_num_active_haplotypes(self):
        return(self.num_active_haplotypes)
    
    
    cdef uint32_t get_num_sites(self, chromosome=None):
        return(self.num_sites)

    
    cdef set_chromosome(self, chromosome):
        """
        Gets alignment index info for this chromosome
        and read the compressed haplotypes metadata.
        """
        if self.current_chromosome != '':
            fclose(self.compressed_tpbwt_file)
        self.current_chromosome = chromosome

        # read in metadata
        md_path = self.path + str(chromosome) + '_md' + self.file_extension
        self.read_haplotype_data(md_path)
        
        # start reading in the compressed haplotypes
        cdef uint32_t num_active_haplotypes
        cdef uint32_t num_all_haplotypes
        cdef uint32_t num_sites
        file_path = self.path + str(chromosome) + self.file_extension
        self.compressed_tpbwt_file = fopen(file_path.encode(), "rb")
        fread(&num_active_haplotypes, sizeof(uint32_t), 1, self.compressed_tpbwt_file)
        self.active_haplotype_indices = <long *>malloc(num_active_haplotypes * sizeof(long))
        fread(&self.active_haplotype_indices[0], sizeof(long), num_active_haplotypes, 
              self.compressed_tpbwt_file)
        fread(&num_all_haplotypes, sizeof(uint32_t), 1, self.compressed_tpbwt_file)
        fread(&num_sites, sizeof(uint32_t), 1, self.compressed_tpbwt_file)
        if num_sites != self.num_sites:
            raise Exception('Error reading compressed TPBWT in ' + file_path
                            + ': wrong number of sites in haplotype alignment.')
        if num_all_haplotypes != self.num_all_haplotypes:
            raise Exception('Error reading compressed TPBWT in ' + file_path
                            + ': wrong number of total haplotypes in alignment.')
        self.num_active_haplotypes = num_active_haplotypes
        self.current_allele_count = 0
        self.haplotype_count = 0


    cdef FILE* get_tpbwt_file(self):
        return(self.compressed_tpbwt_file)
