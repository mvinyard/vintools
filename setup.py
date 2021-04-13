from setuptools import setup
import re
import os
import sys


setup(
    name='vintools',
    version="0.0.2",
    python_requires='>3.8.0',
    author='Michael E. Vinyard - Harvard University - Massachussetts General Hospital - Broad Institute of MIT and Harvard',
    author_email='mvinyard@broadinstitute.org',
    url = 'https://github.com/mvinyard/vintools',
    long_description = open("README.md", encoding="utf-8").read(),
    long_description_content_type = 'text/markdown',
    description="vintools: Michael Vinyard's utilities.",
    packages = [
        'vintools',
	'vintools.math',
        'vintools.machine_learning',
        'vintools.plotting',
        'vintools.utilities',
	'vintools.simulation',
	'vintools.single_cell',
	'vintools.testing_and_development',
    ],
    install_requires=[
        "anndata>=0.7.1",
        "scanpy>=1.4.3",
        "torch>=1.1.0",
        "numpy>=1.19.2",
        "pandas>=1.1.2",
        "pysam>=0.16.0",
        "torchdiffeq>=0.2.1",
    ],
    classifiers = [
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license="MIT"
)

