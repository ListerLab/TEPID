from distutils.core import setup

setup(
    name = 'locaTE',
    version = '1.0',
    description = 'Discover transposable element insertion sites',
    author = 'Tim Stuart',
    author_email = 'timstuart90@gmail.com',
    url = 'https://github.com/timoast/locaTE',
    scripts = ["Scripts/pylocate", "Scripts/locate_map.sh"],
    packages = [''],
)
