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

setup(
    name=name,
    use_scm_version=True,
    install_requires=['Cython', 'pandas', 'numpy'],
    setup_requires=['Cython', 'setuptools-scm'],
    author='23andMe Research',
    author_email='willfreyman@gmail.com',
    description=description,
    license=license,
    long_description=readme,
    long_description_content_type='text/markdown',
    packages=[name],
    include_package_data=True,
    url='https://github.com/23andme/' + name,
    ext_modules=cythonize([Extension('phasedibd.*', ['phasedibd/*.pyx'], include_dirs=[np.get_include()])], language_level="3"),
    include_dirs=[np.get_include()]
)
