#!/usr/bin/env python
from setuptools import setup
import versioneer


setup(
    name = 'TEpy',
    version = versioneer.get_version(),
    description="TEpy: Discover transposable element insertion sites",
    author = 'Tim Stuart',
    author_email = 'timstuart90@gmail.com',
    url = 'https://github.com/timoast/TEpy',
    scripts = ["Scripts/tepy", "Scripts/tepy_map.sh"],
    packages = ['tepy'],
)
