#!/usr/bin/env python
from setuptools import setup
import versioneer
import subprocess


def samtools():
    major, minor = subprocess.check_output(['samtools', '--version']).split()[1].split('.')
    if int(major) >= 1 and int(minor) >= 1:
        return True
    return False


def bedtools():
    p = subprocess.check_output(['bedtools', '--version'])
    major, minor, micro = p.split()[1].split('.')
    major = major[1]
    if int(major) >= 2 and int(minor) >= 24 and int(micro) >= 0:
        return True
    return False


if __name__ == "__main__":
    if not samtools():
        raise Exception("TEpy requires samtools v1.1 or greater")
    if not bedtools():
        raise Exception("TEpy requires bedtools v2.24.0 or greater")


setup(
    name = 'TEpy',
    version = versioneer.get_version(),
    description="TEpy: Discover transposable element insertion sites",
    author = 'Tim Stuart',
    install_requires = [
        'pysam>=0.8.2.1',
        'numpy>=1.9.2',
        'pybedtools>=0.6.9',
    ],
    author_email = 'timstuart90@gmail.com',
    url = 'https://github.com/timoast/TEpy',
    scripts = ["Scripts/tepy_map.sh", "tepy/tepy.py"],
    packages = ['tepy'],
)
