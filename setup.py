from distutils.core import setup
import os

version_file = open(os.path.join(os.path.dirname(__file__), 'VERSION'))
version = version_file.read().strip()

setup(
    name = 'locaTE',
    version = version,
    description = 'Discover transposable element insertion sites',
    author = 'Tim Stuart',
    author_email = 'timstuart90@gmail.com',
    url = 'https://github.com/timoast/locaTE',
    scripts = ["Scripts/pylocate", "Scripts/locate_map.sh"],
    packages = [''],
)
