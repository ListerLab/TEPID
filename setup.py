#!/usr/bin/env python
from setuptools import setup
import versioneer


setup(
    name = 'TEpy',
    version = versioneer.get_version(),
    description="TEpy: Discover transposable element insertion sites",
    author = 'Tim Stuart',
    install_requires = [
        'pysam>=0.8.2.1',
        'numpy>=1.9.2',
        'pybedtools>=0.8.2.1',
    ],
    author_email = 'timstuart90@gmail.com',
    url = 'https://github.com/timoast/TEpy',
    scripts = ["Scripts/tepy_map.sh", "tepy/tepy.py"],
    packages = ['tepy'],
)
