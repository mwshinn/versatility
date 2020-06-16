# Copyright 2019 Max Shinn <maxwell.shinn@yale.edu>
# 
# This file is part of versatility, available under the GNU GPLv3.

from setuptools import setup

with open("README.md", "r") as f:
    long_desc = f.read()

setup(
    name = 'versatility',
    version = '1.0.1',
    description = 'Versatility - find how closely a node in a graph is associated with a community',
    long_description = long_desc,
    long_description_content_type='text/markdown',
    author = 'Max Shinn',
    author_email = 'maxwell.shinn@yale.edu',
    maintainer = 'Max Shinn',
    maintainer_email = 'maxwell.shinn@yale.edu',
    license = 'GPL3',
    python_requires='>=3.5',
    url='https://github.com/mwshinn/versatility',
    py_modules = ['versatility'],
    install_requires = ['numpy', 'scipy', 'matplotlib', 'bctpy', 'networkx'],
    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Topic :: Scientific/Engineering :: Bio-Informatics']
)
