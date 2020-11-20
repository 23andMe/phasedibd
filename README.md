# 23andMe/phasedibd
## templated positional Burrows-Wheeler transform


![build](https://github.com/23andMe/phasedibd/workflows/build/badge.svg)

`phasedibd` is a Python package developed by the 23andMe Ancestry Research team to compute phase aware identity-by-descent (IBD) 
using the templated positional Burrows-Wheeler transform (TPBWT).
See [here](https://www.biorxiv.org/content/10.1101/2020.09.14.296939v1) for our manuscript describing the algorithm, its parameter options, and its performance in detail.

**contact:** Will Freyman <willf@23andme.com>

## Dependencies

`phasedibd` requires the Python packages `Cython`, `numpy`, and `pandas` to be installed.

## To build, install, and run tests:

After cloning the repo you'll need to compile it. First `cd phasedibd` and then
```
make
python setup.py install 
python tests/unit_tests.py
```
Now in Python you can import the module:
```
import phasedibd as ibd
```

## Parameter options:

Adjust sensitivity to error (and the speed of the analysis) with the `template` argument in the `TPBWTAnalysis` constructor. `template` must be a two dimensional python list; by default it is:    
```
[[1, 0, 1, 0],
 [0, 1, 0, 1],
 [1, 1, 0, 0],
 [0, 0, 1, 1],
 [1, 0, 0, 1],
 [0, 1, 1, 0]]
```
This arrangement of templates guarantees all matches will be found as long as no more than two errors are in any four SNP-window.

When `template` is set to `[[1]]` the TPBWT collapses down to Durbin’s original PBWT. This is fast but highly sensitive to error.

When `template` is set to `[[1,0],[0,1]]` the TPBWT will find all matches as long as no more than one error is found within any two SNP window.

When `template` is set to `[[1,0,0],[0,1,0],[0,0,1]]` the TPBWT will find all matches as long as no more than two errors are found within any three SNP window.

Other parameters are set in the `compute_ibd` method in `TPBWTAnalysis`:

| Parameter | Description |
| --- | --- |
| `chromosome` | Default value `‘all’`; valid values `‘1’`, `‘2’`, ..., `‘22’`, `‘X’`. This is useful for setting up computes that parallelize by chromosome. |
| `segments_out_path` | By default the IBD segments are output as a pandas `DataFrame`. Large analyses will use up all available memory and throw a `MemoryError`. By setting this parameter the IBD segments will instead be written to a file. |
| `use_phase_correction` | Default value is `True`. Turns on/off the phase correction heuristic for the entire analysis. |
| `L_m` | Default value is `300`. The minimum number of SNPs a matching subsegment must span to be included in an IBD segment. |
| `L_f` | Default value is `3.0`. The minimum genetic length (cM) of an IBD segment. |
| `missing_site_threshold` | Default value is `10`. The maximum number of consecutive missing SNPs a match will be extended over. |


## IBD output format:

The output is either a pandas `DataFrame` or a `.csv` file with the following columns.
```
chromosome    start    end    start_cm    end_cm    id1    id2    id1_haplotype    id2_haplotype
1             143      2500   2.375       9.8427    243    634    0                1 
1             3470     7037   12.679      26.348    84     591    1                1                             
```
The `start` and `end` columns are the start and end SNP, not the physical (bp) positions. So the number of SNPs a segment spans is `end` - `start`. The columns `id1_haplotype` and `id2_haplotype` must always be 0 or 1.

Any self IBD is included in the output like this:
```
chromosome    start    end    start_cm    end_cm    id1    id2    id1_haplotype    id2_haplotype
1             7467     9940   27.236      50.621    586    586    1                0   
```


For in-sample analyses `id1` < `id2`. For out-of-sample analyses `id1` is always from the first haplotype alignment and `id2` is from the second haplotype alignment.

## To run a basic in-sample IBD analysis on VCF files:

First import the module:
```
import phasedibd as ibd
```
Now create an object to hold the haplotype alignment. This object parses the VCF file;
it expects a VCF file with a diploid GT field. In case of haploid data, the GT field must be
transformed to a pseudo-diploid field (such as 0 -> 0|0).
Additionally, the sites in the VCF file must be sorted by physical position.
There should be one VCF file per chromosome.
```
haplotypes = ibd.VcfHaplotypeAlignment('chr22_sorted.vcf')
```
Next we perform the IBD analysis -- 
instantiate an object of class `TPBWTAnalysis` and call
its method `compute_ibd()`. 
This method has many options (described above) for performing different types of analyses.
The default output is a pandas `DataFrame` with all the IBD segments.
```
tpbwt = ibd.TPBWTAnalysis()
ibd_results = tpbwt.compute_ibd(haplotypes)
```

Optionally you can specify a genetic map:
```
haplotypes = ibd.VcfHaplotypeAlignment('chr22_sorted.vcf', 'chr22_genetic_map.map')
```
If used, it should be a PLINK format genetic map file and contain the same SNPs found in the VCF file. The map file format is:
```
22 . 0.9 16888577
22 . 1.0 16900001
22 . 1.5 17007138
22 . 1.7 19107656
```
If no genetic map is used then the physical positions in the VCF file will be converted to cM assuming 1e6 bp = 1 cM.

## To run a basic out-of-sample IBD analysis:

Runs an IBD analysis between two sets of haplotypes:
```
import phasedibd as ibd
haplotypes_a = ibd.VcfHaplotypeAlignment('a_chr1.vcf', 'chr1_genetic_map.map')
haplotypes_b = ibd.VcfHaplotypeAlignment('b_chr1.vcf', 'chr1_genetic_map.map')
tpbwt = ibd.TPBWTAnalysis()
ibd_results = tpbwt.compute_ibd(haplotypes_a, haplotypes_b)
```
The output is a pandas DataFrame with only the IBD segments shared between the two sets of haplotypes. 
It does not include any IBD shared within the same alignment object.

## Generating TPBWT-compressed haplotypes:

TPBWT-compressed haplotypes are stored in the `.tpbwt` binary file format. 
TPBWT-compressed haplotypes are useful for fast and efficient out-of-sample IBD computes against very large cohort panels.

To TPBWT-compress the haplotypes, set up the haplotypes just like above but now pass them into the `compress_alignment()` method:
```
haplotypes = ibd.VcfHaplotypeAlignment('chr1_sorted.vcf', 'chr1_genetic_map.map')
tpbwt.compress_alignment('compressed_haplotypes/', haplotypes)
```
This will generate files for each chromosome (`1.tpbwt`, `2.tpbwt`, …, `X.tpbwt`) in the directory `compressed_haplotypes/`.

To combine two haplotype alignments into a larger TPBWT-compressed haplotype:
```
tpbwt.compress_alignment('combined_compressed_haplotypes/', 
    haplotypes_1, haplotypes_2)
```

## Generating TPBWT-compressed haplotypes while simultaneously computing IBD:

The `compress_alignment()` method described above is the most memory efficient way to TPBWT-compress haplotypes.
However, it is also possible to simultaneously TPBWT-compress haplotypes and compute IBD.
Set up an IBD compute just like above but use the `compressed_out_path` argument:
```
haplotypes = ibd.VcfHaplotypeAlignment('chr1_sorted.vcf', 'chr1_genetic_map.map')
ibd_segs = tpbwt.compute_ibd(haplotypes, compressed_out_path=’compressed_haplotypes/’)
```
During the IBD compute this will generate files for each chromosome (`1.tpbwt`, `2.tpbwt`, …, `X.tpbwt`) in the directory `compressed_haplotypes/`.

The IBD compute will be slightly slower due to writing the binary `.tpbwt` files.

To combine two haplotype alignments into a larger TPBWT-compressed haplotype:
```
ibd_segs = tpbwt.compute_ibd(haplotypes_1, haplotypes_2, 
    compressed_out_path=’combined_compressed_haplos/’)
```

## Using TPBWT-compressed haplotypes:

Create `CompressedHaplotypeAlignment` objects with the path to the directory
that contains all the TPBWT-compressed files. A simple out-of-sample analysis:
```
haplotypes_1 = ibd.CompressedHaplotypeAlignment(‘compressed_haplotypes_1/’)
haplotypes_2 = ibd.CompressedHaplotypeAlignment(‘compressed_haplotypes_2/’)
tpbwt = ibd.TPBWTAnalysis()
ibd_segs = tpbwt.compute_ibd(haplotypes_1, haplotypes_2)
```
Note that any `CompressedHaplotypeAlignment` objects you want to analyze together must be compressed with the same template and the same set of SNPs.
`CompressedHaplotypeAlignment` objects can be analyzed together with `VcfHaplotypeAlignment` objects
as long as they share the same set of SNPs.

## Other examples:

See the test file for more examples on how to set up different types of analyses: `tests/unit_tests.py`
