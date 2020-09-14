
all:
	pip install cython --user
	python setup.py build_ext --inplace

test:
	python tests/unit_tests.py

clean::
	@find . -name \*~ -exec rm '{}' +
	@find . -name \*.pyc -exec rm '{}' +
	@find . -name __pycache__ -prune -exec rm -vfr '{}' +
	@find phasedibd -name \*.c -prune -exec rm -vfr '{}' +
	@find phasedibd -name \*.so -prune -exec rm -vfr '{}' +
	@rm -rf build bdist cover dist sdist
	@rm -rf .tox .eggs
	@find . \( -name \*.orig -o -name \*.bak -o -name \*.rej \) -exec rm '{}' +
	@rm -rf distribute-* *.egg *.egg-info *.tar.gz cover junit.xml coverage.xml .cache
	@rm -rf compressed_haplotypes_1kgp

.PHONY: clean test
