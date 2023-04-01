from setuptools import setup, find_packages

VERSION = '1.5'
DESCRIPTION = 'Electric field analyses for molecular simulations.'

# Setting up
setup(
    name="tupa",
    version=VERSION,
    author="Marcelo D. Poleto, Justin A Lemkul",
    author_email="Marcelo D. Poleto <mdpoleto@vt.edu>",
    #contributors=CONTRIBUTORS,
    description=DESCRIPTION,
    packages=find_packages(),

    # other arguments...
    #conda_channels=['conda-forge', 'schrodinger'],
    #conda_dependencies=[
    #    'pymol-bundle',
    #],

    install_requires=['MDAnalysis>=2.2.0','numpy', 'configparser', 'argparse', 'warnings'],
    extras_require = {},
    keywords=['python', 'molecular dynamics', 'electric fields', 'MD trajectory analysis'],
    url="https://github.com/mdpoleto/tupa/",
    classifiers=[
        "Programming Language :: Python :: 3.x",
        "License :: OSI Approved :: GNU General Public License v3",
    ],
)
