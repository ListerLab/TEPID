try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup
import os
import subprocess

version_file = os.path.join(os.path.dirname(__file__), 'VERSION')

# update version file if in git repo
try:
    version_git = subprocess.check_output(["git", "describe"]).rstrip()
except:
    version_git = open(version_file).read().strip()
with open(version_file, 'w') as fh:
    fh.write(version_git)


setup(
    name = 'TEpy',
    version = version_git,
    description = 'Discover transposable element insertion sites',
    author = 'Tim Stuart',
    author_email = 'timstuart90@gmail.com',
    url = 'https://github.com/timoast/TEpy',
    scripts = ["Scripts/tepy", "Scripts/tepy_map.sh"],
    packages = [''],
)
