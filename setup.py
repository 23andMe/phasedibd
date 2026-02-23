try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy as np

with open('README.md', 'r') as f:
    readme = f.read()
with open('LICENSE.txt', 'r') as f:
    license = f.read()
name = 'phasedibd'
description = 'Python package for computing the templated positional Burrows-Wheeler transform (TPBWT) used to estimate identity-by-descent.'

extensions = [
    Extension('phasedibd.haplotype_alignment', ['phasedibd/haplotype_alignment.pyx'], include_dirs=[np.get_include()]),
    Extension('phasedibd.compressed_haplotype_alignment', ['phasedibd/compressed_haplotype_alignment.pyx'], include_dirs=[np.get_include()]),
    Extension('phasedibd.vcf_haplotype_alignment', ['phasedibd/vcf_haplotype_alignment.pyx'], include_dirs=[np.get_include()]),
    Extension('phasedibd.templated_pbwt_analysis', ['phasedibd/templated_pbwt_analysis.pyx'], include_dirs=[np.get_include()]),
]

setup(
    name=name,
    use_scm_version=True,
    install_requires=['Cython', 'pandas', 'numpy'],
    author='23andMe Research',
    author_email='willfreyman@gmail.com',
    description=description,
    license=license,
    long_description=readme,
    long_description_content_type='text/markdown',
    packages=[name],
    include_package_data=True,
    url='https://github.com/23andme/' + name,
    ext_modules=cythonize(extensions, include_path=['phasedibd']),
    include_dirs=[np.get_include()]
)
