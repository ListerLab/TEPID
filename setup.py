#!/usr/bin/env python

from setuptools import setup
import versioneer
import subprocess


def samtools():
    v = subprocess.check_output(['samtools', '--version']).split()[1].split('.')
    major = int(v[0])
    minor = int(v[1])
    if major >= 1 and minor >= 4:
        return True
    return False


def bedtools():
    v = subprocess.check_output(['bedtools', '--version']).split()[1].split('.')
    major = v[0].strip('v')
    if int(major) >= 2 and int(v[1]) >= 24:
        return True
    return False


if __name__ == "__main__":
    if not samtools():
        raise Exception("TEPID requires samtools >= v1.1 and < v1.3")
    if not bedtools():
        raise Exception("TEPID requires bedtools v2.25.0 or greater")


setup(
    name = 'TEPID',
    version = versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="TEPID: transposable element polymorphism identification",
    author = 'Tim Stuart',
    install_requires = [
        'pysam<0.9, >0.8',
        'numpy>=1.9.2',
        'pybedtools>=0.6.9',
        'pandas',
        'nose'
    ],
    author_email = 'timstuart90@gmail.com',
    url = 'https://github.com/ListerLab/TEPID',
    scripts = ["Scripts/tepid-map", "Scripts/tepid-map-se", "Scripts/tepid-discover", "Scripts/tepid-refine"],
    packages = ['tepid'],
    test_suite="nose.collector"
)
